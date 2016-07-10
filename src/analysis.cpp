#include "../include/analysis.h"
#include "../include/paramdata.h"
#include "../include/nsga2.h"

// コンストラクタ
Analysis::Analysis(double arg_l, double arg_b, double arg_a, double arg_m1, double arg_m2)
{
	this->_lambda	= arg_l;
	this->_beta2	= arg_b;
	this->_alpha	= arg_a;
	this->mu1	= arg_m1;
	this->mu2	= arg_m2;
}

// デストラクタ
Analysis::~Analysis() {}

void Parameter::setParameter(std::vector<double> &a_, std::vector<double> &mu1_, std::vector<double> &mu2_, std::vector<double> &sigma1_, std::vector<double> &sigma2_, std::vector<double> &kappa_)
{
	std::copy(a_.begin(), a_.end(), back_inserter(_a));
	std::copy(mu1_.begin(), mu1_.end(), back_inserter(_mu1));
	std::copy(mu2_.begin(), mu2_.end(), back_inserter(_mu2));
	std::copy(sigma1_.begin(), sigma1_.end(), back_inserter(_sigma1));
	std::copy(sigma2_.begin(), sigma2_.end(), back_inserter(_sigma2));
	std::copy(kappa_.begin(), kappa_.end(), back_inserter(_kappa));
}

void Parameter::allocParameter()
{
	_a.reserve(NUM_GAUSS);
	_mu1.reserve(NUM_GAUSS);
	_mu2.reserve(NUM_GAUSS);
	_sigma1.reserve(NUM_GAUSS);
	_sigma2.reserve(NUM_GAUSS);
	_kappa.reserve(NUM_GAUSS);
}

void Parameter::freeParameter()
{
	std::vector<double>().swap(_a);
	std::vector<double>().swap(_mu1);
	std::vector<double>().swap(_mu2);
	std::vector<double>().swap(_sigma1);
	std::vector<double>().swap(_sigma2);
	std::vector<double>().swap(_kappa);
}

std::vector<double> Parameter::getParameter(std::string prm)
{
	if (prm == "a")
		return _a;
	else if (prm == "mu1")
		return _mu1;
	else if (prm == "mu2")
		return _mu2;
	else if (prm == "sigma1")
		return _sigma1;
	else if (prm == "sigma2")
		return _sigma2;
	else if (prm == "kappa")
		return _kappa;
	else {
		std::cout << "Error: Given parameter doesn't exist." << std::endl;
		return std::vector<double>();
	}
}

/**
 * @fn NSGA2でモーメント方程式を解く
 * @param Parameter* prm モーメント方程式から求めたパラメータ
 */
int Analysis::GeneticAlgorithm(std::vector<Parameter*> &prm)
{
	std::cout << "Get analysis solution using NSGA2.\n" << std::endl;

	unsigned int i;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*this->_beta2);

	// 入力に関するモーメント
	double dF[6];
	dF[0] = 0;
	dF[1] = this->_alpha*S0 + this->_lambda*(1. - this->_alpha)*this->_beta2;
	dF[2] = 0;
	dF[3] = this->_lambda*pow((1. - this->_alpha), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = this->_lambda*pow((1. - this->_alpha), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));

	// 入力や系のパラメータ
	ParamData* setData = new ParamData(NUM_OF_MOMENTEQ, NUM_OF_PARAM, ZETA, EPSILON, dF);

	// nsga2
	NSGA2 *n2	= new NSGA2(10, 20, true, 300);
	n2->run(setData);

		/********** 計算結果 **********/
	prm.reserve(n2->getPrmValue().size());
	for (i = 0; i < n2->getPrmValue().size(); ++i) {
		std::vector<double> prm_a(NUM_GAUSS), prm_mu1(NUM_GAUSS), prm_mu2(NUM_GAUSS), prm_sigma1(NUM_GAUSS), prm_sigma2(NUM_GAUSS), prm_kappa(NUM_GAUSS);
		// 重み
		prm_a[1] = n2->getPrmValue()[i][0];		// a2
		prm_a[0] = (1. - prm_a[1]) / 2.;		// a1
		prm_a[2] = prm_a[0];					// a3
		// 変位
		prm_mu1[0] = n2->getPrmValue()[i][1];	// mu11
		prm_mu1[1] = 0.;						// mu21
		prm_mu1[2] = -1.*prm_mu1[0];			// mu31
		prm_sigma1[0] = n2->getPrmValue()[i][3];	// sigma11
		prm_sigma1[1] = n2->getPrmValue()[i][5];	// sigma21
		prm_sigma1[2] = prm_sigma1[0];				// sigma31
		// 速度
		prm_mu2[0] = n2->getPrmValue()[i][2];	// mu12
		prm_mu2[1] = 0.;						// mu22
		prm_mu2[2] = -1.*prm_mu2[0];			// mu32
		prm_sigma2[0] = n2->getPrmValue()[i][4];	// sigma12
		prm_sigma2[1] = n2->getPrmValue()[i][6];	// sigma22
		prm_sigma2[2] = prm_sigma2[0];				// sigma32
		// 共分散
		prm_kappa[0] = n2->getPrmValue()[i][7];		// kappa123
		prm_kappa[1] = n2->getPrmValue()[i][8];		// kappa223
		prm_kappa[2] = n2->getPrmValue()[i][9];		// kappa323
	
		Parameter* p	= new Parameter();
		p->allocParameter();
		p->setParameter(prm_a, prm_mu1, prm_mu2, prm_sigma1, prm_sigma2, prm_kappa);
		prm.push_back(p);
	}
	/******************************/

	delete setData;
	delete n2;

	return EXIT_SUCCESS;
}

/**
 * @fn 最小二乗法を解いて解析解を求める
 * @param Parameter* prm モーメント方程式から求めたパラメータ
 */
std::string Analysis::leastSquareMethod(Parameter* prm)
{
	std::cout << "Get analysis solution using least squeare method.\n" << std::endl;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*this->_beta2);

	// 入力に関するモーメント
	double dF[6];
	dF[0] = 0;
	dF[1] = this->_alpha*S0 + this->_lambda*(1. - this->_alpha)*this->_beta2;
	dF[2] = 0;
	dF[3] = this->_lambda*pow((1. - this->_alpha), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = this->_lambda*pow((1. - this->_alpha), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));

	// 等価線形化法で初期値を計算
	double sigma_x, sigma_y, rho_xy;
	this->_culcInitValue(&sigma_x, &sigma_y, &rho_xy);

	/*  初期値         {a,   μ1,           μ2,          σ11,     σ12,     σ21,     σ22,     k1,                    k2, k3}*/
	double x_init[NUM_OF_PARAM] = { 0.5, sigma_x + this->mu1, sigma_y + this->mu2, sigma_x, sigma_y, sigma_x, sigma_y, rho_xy*sigma_x*sigma_y, 0., 0. };
	
	// 最小二乗法で使うパラメータ
	ParamData* setData = new ParamData(NUM_OF_MOMENTEQ, NUM_OF_PARAM, ZETA, EPSILON, dF);

	// 最小二乗法を解くための関数をセット
	gsl_multifit_function_fdf f;
	f.f		= &MomentEq::expb_f;
	f.df	= &MomentEq::expb_df;
	f.fdf	= &MomentEq::expb_fdf;
	f.n		= NUM_OF_MOMENTEQ;
	f.p		= NUM_OF_PARAM;
	f.params = setData;

	// 最小二乗法のソルバーをセット
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, NUM_OF_MOMENTEQ, NUM_OF_PARAM);

	// 計算開始
	gsl_vector_view x;
	x = gsl_vector_view_array(x_init, NUM_OF_PARAM);
	gsl_multifit_fdfsolver_set(s, &f, &x.vector);
	size_t iter = 0;
	int status;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if (status) break;
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-10, 1e-10);
	} while (status == GSL_CONTINUE && iter < 10000);


	/********** 計算結果 **********/
	std::vector<double> prm_a(NUM_GAUSS), prm_mu1(NUM_GAUSS), prm_mu2(NUM_GAUSS), prm_sigma1(NUM_GAUSS), prm_sigma2(NUM_GAUSS), prm_kappa(NUM_GAUSS);
	// 重み
	prm_a[1] = gsl_vector_get(s->x, 0);		// a2
	prm_a[0] = (1. - prm_a[1]) / 2.;		// a1
	prm_a[2] = prm_a[0];					// a3
	// 変位
	prm_mu1[0] = gsl_vector_get(s->x, 1);	// mu11
	prm_mu1[1] = 0.;						// mu21
	prm_mu1[2] = -1.*prm_mu1[0];			// mu31
	prm_sigma1[0] = gsl_vector_get(s->x, 3);	// sigma11
	prm_sigma1[1] = gsl_vector_get(s->x, 5);	// sigma21
	prm_sigma1[2] = prm_sigma1[0];				// sigma31
	// 速度
	prm_mu2[0] = gsl_vector_get(s->x, 2);	// mu12
	prm_mu2[1] = 0.;						// mu22
	prm_mu2[2] = -1.*prm_mu2[0];			// mu32
	prm_sigma2[0] = gsl_vector_get(s->x, 4);	// sigma12
	prm_sigma2[1] = gsl_vector_get(s->x, 6);	// sigma22
	prm_sigma2[2] = prm_sigma2[0];				// sigma32
	// 共分散
	prm_kappa[0] = gsl_vector_get(s->x, 7);		// kappa123
	prm_kappa[1] = gsl_vector_get(s->x, 8);		// kappa223
	prm_kappa[2] = gsl_vector_get(s->x, 9);		// kappa323

	prm->setParameter(prm_a, prm_mu1, prm_mu2, prm_sigma1, prm_sigma2, prm_kappa);
	/******************************/

	// TODO:何をやっているのか見直す
	gsl_matrix *J = gsl_matrix_alloc(NUM_OF_MOMENTEQ, NUM_OF_PARAM);
	gsl_multifit_fdfsolver_jac(s, J);
	gsl_matrix *covar = gsl_matrix_alloc(NUM_OF_PARAM, NUM_OF_PARAM);
	gsl_multifit_covar(J, 0.0, covar);

#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	double chi = gsl_blas_dnrm2(s->f);
	double dof = NUM_OF_MOMENTEQ - NUM_OF_PARAM;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));

	std::cout << "roop(iter): " << iter << std::endl;
	std::cout << "/*********************solution*********************/" << std::endl;
	std::cout << "a	      = " << prm->getParameter("a")[1] << " +/- " << c*ERR(0) << std::endl;
	std::cout << "mu1     = " << prm->getParameter("mu1")[0] << " +/- " << c*ERR(1) << std::endl;
	std::cout << "mu2     = " << prm->getParameter("mu2")[0] << " +/- " << c*ERR(2) << std::endl;
	std::cout << "sigma11 = " << prm->getParameter("sigma1")[0] << " +/- " << c*ERR(3) << std::endl;
	std::cout << "sigma11 = " << prm->getParameter("sigma2")[0] << " +/- " << c*ERR(4) << std::endl;
	std::cout << "sigma21 = " << prm->getParameter("sigma1")[1] << " +/- " << c*ERR(5) << std::endl;
	std::cout << "sigma22 = " << prm->getParameter("sigma2")[1] << " +/- " << c*ERR(6) << std::endl;
	std::cout << "k1      = " << prm->getParameter("kappa")[0] << " +/- " << c*ERR(7) << std::endl;
	std::cout << "k2      = " << prm->getParameter("kappa")[1] << " +/- " << c*ERR(8) << std::endl;
	std::cout << "k3      = " << prm->getParameter("kappa")[2] << " +/- " << c*ERR(9)<< std::endl;
	std::cout << "/**************************************************/" << std::endl;
	std::cout << "status  = " << gsl_strerror(status) << "\n" << std::endl;

	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free (J);

	delete setData;

	return gsl_strerror(status);
}

/**
 * @fn 変位のPDFを求める
 * @param Parameter* prm モーメント方程式から求めたパラメータ
 * @param vector &x X軸の情報を保存（領域確保済み）
 * @param vector &y Y軸の情報を保存（領域確保済み）
 * @param int xmin X軸の最小値
 */
void Analysis::createDispPdf(Parameter* prm, std::vector<double> &x, std::vector<double> &y, int xmin)
{
	std::cout << "Culculate the displacement pdf.\n" << std::endl;

	unsigned int i;
	double integration = 0.;

	for (i = 0; i*0.01 <= abs(xmin)*2; ++i) {
		x[i]	= xmin + i*0.01;
		y[i]	= this->_createGaussianPdf(prm->getParameter("a"), prm->getParameter("mu1"), prm->getParameter("sigma1"), xmin + i*0.01);
		integration += 0.01*this->_createGaussianPdf(prm->getParameter("a"), prm->getParameter("mu1"), prm->getParameter("sigma1"), xmin + i*0.01);
	}

	std::cout << "y1 integration = " << integration << ".\n" << std::endl;
}

// 変位の尖り度を求める
void Analysis::culcDispVarience()
{
	std::cout << "Creating a file of the displacement varience(.dat).\n" << std::endl;

	FILE *y1_var, *gsax1pdf;
	y1_var = fopen("anl_y1_var.dat", "w");
	gsax1pdf = fopen("gsay1pdf.dat", "r");

	double var_y1 = 0., mo4_y1 = 0.;
	double row1, row2;
	int line = 0, ret1;
	while ((ret1 = fscanf(gsax1pdf, "%lf %lf", &row1, &row2)) != EOF)
	{
		if ((line % 10) == 0)
		{
			var_y1 += pow(row1, 2.)*row2;
			mo4_y1 += pow(row1, 4.)*row2;
		}
		line++;
	}

	fprintf(y1_var, "%lf", mo4_y1 / pow(var_y1, 2.));

	fclose(gsax1pdf);
	fclose(y1_var);
}

// 閾値通過率の分布を求める
void Analysis::createLevelCrossing(Parameter *prm)
{
	std::cout << "Creating a file of the level crossing(.dat).\n" << std::endl;

	unsigned int xi;

	FILE *x1_prob_pass;
	x1_prob_pass = fopen("x1_prob_pass.dat", "w");

	double pp_xi, prob_pass = 0.;

	for (xi = 0; xi*0.01 < 8; ++xi)
	{
		prob_pass = 0.;	// 元に戻す
		pp_xi = (double)xi*0.01;

		prob_pass = this->_culcLevelCrossing(pp_xi, prm->getParameter("a"), prm->getParameter("mu1"), prm->getParameter("mu2"), prm->getParameter("sigma1"), prm->getParameter("sigma2"), prm->getParameter("kappa"));

		fprintf(x1_prob_pass, "%lf %lf\n", pp_xi, prob_pass);
	}

	fclose(x1_prob_pass);
}

/**
 * @fn 速度のPDFを求める
 * @param Parameter* prm モーメント方程式から求めたパラメータ
 * @param vector &x X軸の情報を保存（領域確保済み）
 * @param vector &y Y軸の情報を保存（領域確保済み）
 * @param int xmin X軸の最小値
 */
void Analysis::createVelPdf(Parameter *prm, std::vector<double> &x, std::vector<double> &y, int xmin)
{
	std::cout << "Creating a file of the velocity varience(.dat).\n" << std::endl;

	unsigned int i;
	double integration = 0.;

	for (i = 0; i*0.01 <= 30; ++i)
	{
		x[i]	= xmin + i*0.01;
		y[i]	= this->_createGaussianPdf(prm->getParameter("a"), prm->getParameter("mu2"), prm->getParameter("sigma2"), xmin);
		integration += 0.01*this->_createGaussianPdf(prm->getParameter("a"), prm->getParameter("mu2"), prm->getParameter("sigma2"), xmin);
	}
	printf("y integration = %lf \n",integration);
}

/**
 * 等価線形化法により初期値を計算
 * 線形の連立方程式Ax=Bを解く
 */
void Analysis::_culcInitValue(double *sigma_x, double *sigma_y, double *rho_xy)
{
	// カウント変数
	unsigned int tmp;
	
	int s, num = 3;
	double a[num*num], b[num];
	double ke, Exxold, Exyold, Eyyold, Exx=0.5, Exy=0., Eyy=0., err = 1000.;	// ループ計算時の解の保存用，x:変位，y:速度

	/**
		応答のガウス性を仮定し，等価線形化した系のモーメント方程式を解く
			2*Exy = 0
			-2*ke*Exx - 2*ZETA*Exy + Eyy = 0
			-2*Exy -4*ZETA*Eyy + (alpha*2*PI*S0 + (1-alpha)*lambda*beta2)/dt = 0
	 */
	for (tmp = 0; err >10e-6; tmp++)
	{
		// ひとつ前の情報を保存
		Exxold = Exx;
		Exyold = Exy;
		Eyyold = Eyy;

		ke = 1.+3.*EPSILON*Exx;	// 等価線形係数 ke = 1+3εE[X^2]

		a[0*num+0]	= 0.;	a[0*num+1]	= 2.;		a[0*num+2]	= 0.;
		a[1*num+0]	= -ke;	a[1*num+1]	= 2.*ZETA;	a[1*num+2]	= 1.;
		a[2*num+0]	= 0.;	a[2*num+1]	= 2.*ke;	a[2*num+2]	= 4.*ZETA;
		b[0]	= 0.;
		b[1]	= 0.;
		b[2]	= this->_alpha*2.*PI*S0 + (1.-this->_alpha)*this->_lambda*this->_beta2;
		
		// LU分解の方法でモーメント方程式を解く
		gsl_matrix_view m	= gsl_matrix_view_array(a, num, num);
		gsl_vector_view c	= gsl_vector_view_array(b, num);
		gsl_vector *x		= gsl_vector_alloc(num);
		gsl_permutation *px	= gsl_permutation_alloc(num);
		gsl_linalg_LU_decomp(&m.matrix, px, &s);
		gsl_linalg_LU_solve(&m.matrix, px, &c.vector, x);
		Exx = gsl_vector_get(x,0);
		Exy = gsl_vector_get(x,1);
		Eyy = gsl_vector_get(x,2);
		gsl_vector_free(x);
		gsl_permutation_free(px);

		if (tmp == 0) continue;
		err = pow(Exx - Exxold, 2.) + pow(Exy - Exyold, 2.) + pow(Eyy - Eyyold, 2.);	// 収束条件に使う誤差，前ループとの差の二乗和
	}
	
	*sigma_x = sqrt(Exx);
	*sigma_y = sqrt(Exy);
	*rho_xy  = Exy/sqrt(Exx*Eyy);

	std::cout << "roop(_culcInitValue):" <<  tmp << std::endl;
	std::cout << "sigma_x=" << sqrt(Exx) << ", sigma_y=" << sqrt(Eyy) << ", rho_xy=" << Exy/sqrt(Exx*Eyy) << "\n" << std::endl;
	return;
}

/**
 * @fn ファイルに出力する
 * @param string name ファイル名
 * @param vector &x X軸情報
 * @param vector &y Y軸情報
 */
void Analysis::outputIntoFile(const std::string name, const std::vector<double> &x, const std::vector<double> &y)
{
	std::cout << "Creating a file.\n" << std::endl;

	// 値の数をチェック
	if (x.size() != y.size()) {
		std::cout << "Error: Do not match vector size." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ofstream ofs(name);
	// FILE *f;
	// if ((f = fopen(name, "w")) == NULL) {
	// 	std::cout << "Error: Failed to open file." << std::endl;
	// 	exit(EXIT_FAILURE);
	// }

	unsigned int i;
	for (i = 0; i < x.size(); ++i) {
		ofs << x[i] << " " << y[i] << std::endl;
		// fprintf(f, "%lf %lf\n", x[i], y[i]);
	}
	// fclose(f);
}

/**
 * 1変数のガウス分布を足し合わせる
 * 
 * @param double a[] 重み
 * @param double mu[] 平均
 * @param double sigma[] 分散
 * @param double x 変数
 * @return pdf 確率密度関数
 */
double Analysis::_createGaussianPdf(const std::vector<double> &a, const std::vector<double> &mu, const std::vector<double> &sigma, double x) 
{
	double pdf = 0.;
	unsigned int i;

	for (i = 0; i < NUM_GAUSS; ++i)
		pdf += a[i]*(1./sqrt(2.*PI)/sigma[i])*exp(-(x-mu[i])*(x-mu[i])/2./pow(sigma[i],2.));

	return pdf;
}

/**
 * 閾値通過率を計算
 * @param double pp_xi 閾値
 * @param double a[] 重み
 * @param double mu1[] 変位の平均
 * @param double mu2[] 速度の平均
 * @param double sigma1[] 変位の分散
 * @param double sigma2[] 速度の分散
 * @param double kappa[]
 */
double Analysis::_culcLevelCrossing(double pp_xi, const std::vector<double> &a, const std::vector<double> &mu1, const std::vector<double> &mu2, const std::vector<double> &sigma1, const std::vector<double> &sigma2, const std::vector<double> &kappa)
{
	double prob_pass = 0.;
	double pp_c = 0., pp_g = 0., pp_sigma = 0.;
	int tmp;

	for (tmp = 0; tmp < NUM_GAUSS; ++tmp)
	{
		pp_c		= kappa[tmp]/sigma1[tmp]/sigma2[tmp];
		pp_g		= mu2[tmp] + pp_c*sigma2[tmp]*(pp_xi - mu1[tmp])/sigma1[tmp];
		pp_sigma	= sigma2[tmp]*sqrt(1. - pow(pp_c,2.));
		// 閾値通過率
		prob_pass	+= a[tmp]*exp(-pow((pp_xi - sigma1[tmp]),2.)/2./pow(sigma1[tmp],2.)) /2./PI/sigma1[tmp]/sigma2[tmp]/sqrt(1. - pow(pp_c,2.))* (pow(pp_sigma,2.)*exp(-pow(pp_g,2.)/2./pow(pp_sigma,2.))
						+ sqrt(PI/2.)*pp_g*pp_sigma*(1. + erf(pp_g/sqrt(2.)/pp_sigma)));
	}

	return prob_pass;
}
