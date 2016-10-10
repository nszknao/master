#include "../include/analysis.h"
#include "../include/nsga2.h"
#include "../include/nsga3.h"
#include "../include/common.h"

Analysis::Analysis(double arg_l, double arg_b, double arg_a, double arg_m1, double arg_m2)
{
	this->_lambda	= arg_l;
	this->_beta2	= arg_b;
	this->_alpha	= arg_a;
	this->_mu1	= arg_m1;
	this->_mu2	= arg_m2;
}

// デストラクタ
Analysis::~Analysis() {}

/**
 * @fn NSGA2でモーメント方程式を解く
 * @param vector<double> &pValue 計算結果のパラメータ値を保存
 * @param vector<double> &oValue 計算結果の目的関数値を保存
 */
int Analysis::GeneticAlgorithm(std::vector<GAIndividual> &pops)
{
	std::cout << "Get analysis solution using NSGA2.\n" << std::endl;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*this->_beta2);

	// 入力に関するモーメント
	std::vector<double> dF(6);
	dF[0] = 0;
	dF[1] = this->_alpha*S0 + this->_lambda*(1. - this->_alpha)*this->_beta2;
	dF[2] = 0;
	dF[3] = this->_lambda*pow((1. - this->_alpha), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = this->_lambda*pow((1. - this->_alpha), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));

	// nsga2
	// NSGA2 *n2	= new NSGA2(250, 200);
	// n2->run(&dF[0]);
	// pops	= n2->getFinalPops();
	// delete n2;

	// nsga3
	NSGA3 *n3	= new NSGA3();
	n3->run();
	pops	= n3->getFinalPops();
	delete n3;

	return EXIT_SUCCESS;
}

/**
 * @fn 最小二乗法を解いて解析解を求める
 * @param Parameter* prm モーメント方程式から求めたパラメータ
 */
std::string Analysis::leastSquareMethod(std::vector<double> &pValue)
{
	std::cout << "Get analysis solution using least squeare method.\n" << std::endl;

	unsigned int i;
	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*this->_beta2);

	// 入力に関するモーメント
	std::vector<double> dF(6);
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
	double x_init[NUM_OF_PARAM] = { 0.5, sigma_x + this->_mu1, sigma_y + this->_mu2, sigma_x, sigma_y, sigma_x, sigma_y, rho_xy*sigma_x*sigma_y, 0., 0. };
	
	// 最小二乗法を解くための関数をセット
	gsl_multifit_function_fdf f;
	f.f		= &MomentEq::expb_f;
	f.df	= &MomentEq::expb_df;
	f.fdf	= &MomentEq::expb_fdf;
	f.n		= NUM_OF_MOMENTEQ;
	f.p		= NUM_OF_PARAM;
	f.params = &dF[0];

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

	// 計算結果
	pValue.resize(NUM_OF_PARAM);
	for (i = 0; i < NUM_OF_PARAM; ++i) {
		pValue[i]	= gsl_vector_get(s->x, i);
	}

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
	std::cout << "a	      = " << pValue[0] << " +/- " << c*ERR(0) << std::endl;
	std::cout << "mu1     = " << pValue[1] << " +/- " << c*ERR(1) << std::endl;
	std::cout << "mu2     = " << pValue[2] << " +/- " << c*ERR(2) << std::endl;
	std::cout << "sigma11 = " << pValue[3] << " +/- " << c*ERR(3) << std::endl;
	std::cout << "sigma12 = " << pValue[4] << " +/- " << c*ERR(4) << std::endl;
	std::cout << "sigma21 = " << pValue[5] << " +/- " << c*ERR(5) << std::endl;
	std::cout << "sigma22 = " << pValue[6] << " +/- " << c*ERR(6) << std::endl;
	std::cout << "k1      = " << pValue[7] << " +/- " << c*ERR(7) << std::endl;
	std::cout << "k2      = " << pValue[8] << " +/- " << c*ERR(8) << std::endl;
	std::cout << "k3      = " << pValue[9] << " +/- " << c*ERR(9)<< std::endl;
	std::cout << "/**************************************************/" << std::endl;
	std::cout << "status  = " << gsl_strerror(status) << "\n" << std::endl;

	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free (J);

	return gsl_strerror(status);
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
 * @fn ファイルに個体情報を出力する
 * --- 出力形式 ---
 * # E[y1^2] E[y2y1] E[y2^2] ... E[y1^3y2^5] 目的関数1の値 目的関数2の値 ...  目的関数15の値
 * 0.0001 1.4356
 * ...
 * --- 出力形式 ---
 * @param string name ファイル名
 * @param GAIndividual &ind 個体
 * @param vector &x X軸情報
 * @param vector &y Y軸情報
 */
void Analysis::outputPopsIntoFile(const std::string name, const GAIndividual &ind, const std::vector<double> &x, const std::vector<double> &y)
{
	// 値の数をチェック
	if (x.size() != y.size()) {
		std::cout << "Error: Do not match vector size." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ofstream ofs(name);
	unsigned int i;
	ofs << "#";
	for (i = 0; i < ind.mValue.size(); ++i) {
		ofs << " " << ind.mValue[i] << std::flush;
	}
	for (i = 0; i < ind.oValue.size(); ++i) {
		ofs << " " << ind.oValue[i] << std::flush;
	}
	ofs << std::endl;

	for (i = 0; i < x.size(); ++i) {
		ofs << x[i] << " " << y[i] << std::endl;
	}
}

/**
 * @fn ファイルに全個体情報を出力する
 * --- 出力形式 ---
 * # E[y1^2] E[y2y1] E[y2^2] ... E[y1^3y2^5] 目的関数1の値 目的関数2の値 ...  目的関数15の値
 * 0.0001 1.4356
 * ...
 * --- 出力形式 ---
 * @param string name ファイル名
 * @param GAIndividual &ind 個体
 * @param vector &x X軸情報
 * @param vector &y Y軸情報
 */
void Analysis::outputAllPopsIntoFile(const std::string name, const std::vector<GAIndividual> &pops)
{
	unsigned int i, ii;

	std::ofstream ofs(name);
	for (i = 0; i < pops.size(); ++i) {
		for (ii = 0; ii < pops[i].mValue.size(); ++ii) {
			ofs << pops[i].mValue[ii] << " " << std::flush;
		}
		for (ii = 0; ii < pops[i].oValue.size(); ++ii) {
			ofs << pops[i].oValue[ii] << " " << std::flush;
		}
		for (ii = 0; ii < pops[i].pValue.size(); ++ii) {
			if (ii == pops[i].pValue.size()-1) {
				ofs << pops[i].pValue[ii] << std::flush;				
			} else {
				ofs << pops[i].pValue[ii] << " " << std::flush;				
			}
		}
		ofs << std::endl;
	}
}

/**
 * main
 */
int main(int argc, char *argv[])
{
	std::string filename;

	char *ends;
	double lambda	= strtod(argv[1],&ends);
	double beta2	= strtod(argv[2],&ends);
	double alpha	= strtod(argv[3],&ends);
	double mu1		= strtod(argv[4],&ends);


	std::cout << "--------------------\n" << std::endl;
	std::cout << "analysis.cpp started.\n" << std::endl;
	
	Analysis *ana	= new Analysis(lambda, beta2, alpha, mu1, mu1);

	/* 最小二乗法で解く */
	// std::vector<double> pValue;
	// std::string result	= ana->leastSquareMethod(pValue);
	// if (result == "success") {
	// 	Parameter* prm	= new Parameter();
	// 	prm->allocParameter();
	// 	ana->getDetailParameterFromSimpleNotation(prm, pValue);

	// 	Common::outputIntoFile((char*)"gsay1pdf.dat", dispX, dispY);

	// 	delete prm;
	// }

	/* GAで解く */
	std::vector<GAIndividual> pops;
	ana->GeneticAlgorithm(pops);
	filename	= "ana_gsay1pdf.dat";
	ana->outputAllPopsIntoFile(filename, pops);

	delete ana;
	
	std::cout << "analysis.cpp has done.\n" << std::endl;
	std::cout << "--------------------\n" << std::endl;
	return 0;
}