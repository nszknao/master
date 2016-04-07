#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include "expfit.h"
#include "ParamData.h"

#define N 15
#define P 10
#define PI M_PI
#define iN 3		// 初期値を求めるときの方程式の数
#define NUM_GAUSS 3	// 足しあわせるガウス分布の数

/************ 系の係数・入力条件（不変）*******************/
#define S0 1./(2.*PI)
#define ZETA 0.2
#define EPSILON 0.3
#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布


class Analysis
{
private:
	// 混合ガウス分布のパラメータ
	double prm_a[NUM_GAUSS], prm_mu1[NUM_GAUSS], prm_mu2[NUM_GAUSS], prm_sigma1[NUM_GAUSS], prm_sigma2[NUM_GAUSS], prm_kappa[NUM_GAUSS];
	double _createGaussianPdf(double a[], double mu[], double sigma[], double x);
	double _culcLevelCrossing(double pp_xi, double a[], double mu1[], double mu2[], double sigma1[], double sigma2[], double kappa[]);
	void _culcInitValue(double lambda, double beta2, double ggd_kappa, double alpha, double *sigma_x, double *sigma_y, double *rho_xy);

public:
	double lambda, beta2, alpha, mu1, mu2;
	std::string leastSquareMethod();
	void createDispPdf();
	void createVelPdf();
	void culcDispVarience();
	void culcVelVarience();
	void createLevelCrossing();
};

// 最小二乗法を解いて解析解を求める
std::string Analysis::leastSquareMethod()
{
	std::cout << "Get analysis solution using least squeare method.\n" << std::endl;

	// カウント変数
	size_t tmp;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*beta2);

	// 入力に関するモーメント
	double dF[6];
	dF[0] = 0;
	dF[1] = this->alpha*S0 + this->lambda*(1. - this->alpha)*(pow(ggd_a, 2.)*gsl_sf_gamma(3. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[2] = 0;
	dF[3] = this->lambda*pow((1. - alpha), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = this->lambda*pow((1. - alpha), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));

	/* 近似対象となる観測データを生成 */
	double y[N];
	for (tmp = 0; tmp < N; tmp++) y[tmp] = 0.;

	// 等価線形化法で初期値を計算
	double sigma_x, sigma_y, rho_xy;
	this->_culcInitValue(this->lambda, this->beta2, GGD_KAPPA, this->alpha, &sigma_x, &sigma_y, &rho_xy);

	/*  初期値         {a,   μ1,           μ2,          σ11,     σ12,     σ21,     σ22,     k1,                    k2, k3}*/
	double x_init[P] = { 0.5, sigma_x + this->mu1, sigma_y + this->mu2, sigma_x, sigma_y, sigma_x, sigma_y, rho_xy*sigma_x*sigma_y, 0., 0. };
	
	// 最小二乗法で使うパラメータ
	ParamData* setData = new ParamData(N, P, y, ZETA, EPSILON, dF);

	// 最小二乗法を解くための関数をセット
	gsl_multifit_function_fdf f;
	f.f		= &MomentEq::expb_f;
	f.df	= &MomentEq::expb_df;
	f.fdf	= &MomentEq::expb_fdf;
	f.n		= N;
	f.p		= P;
	f.params = setData;

	// 最小二乗法のソルバーをセット
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, N, P);

	// 計算開始
	gsl_vector_view x;
	x = gsl_vector_view_array(x_init, P);
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
	// 重み
	this->prm_a[1] = gsl_vector_get(s->x, 0);		// a2
	this->prm_a[0] = (1. - this->prm_a[1]) / 2.;		// a1
	this->prm_a[2] = this->prm_a[0];					// a3
	// 変位
	this->prm_mu1[0] = gsl_vector_get(s->x, 1);	// mu11
	this->prm_mu1[1] = 0.;						// mu21
	this->prm_mu1[2] = -1.*this->prm_mu1[0];			// mu31
	this->prm_sigma1[0] = gsl_vector_get(s->x, 3);	// sigma11
	this->prm_sigma1[1] = gsl_vector_get(s->x, 5);	// sigma21
	this->prm_sigma1[2] = this->prm_sigma1[0];				// sigma31
	// 速度
	this->prm_mu2[0] = gsl_vector_get(s->x, 2);	// mu12
	this->prm_mu2[1] = 0.;						// mu22
	this->prm_mu2[2] = -1.*this->prm_mu2[0];			// mu32
	this->prm_sigma2[0] = gsl_vector_get(s->x, 4);	// sigma12
	this->prm_sigma2[1] = gsl_vector_get(s->x, 6);	// sigma22
	this->prm_sigma2[2] = this->prm_sigma2[0];				// sigma32
	// 共分散
	this->prm_kappa[0] = gsl_vector_get(s->x, 7);		// kappa123
	this->prm_kappa[1] = gsl_vector_get(s->x, 8);		// kappa223
	this->prm_kappa[2] = gsl_vector_get(s->x, 9);		// kappa323
	/******************************/

	// TODO:何をやっているのか見直す
	gsl_matrix *covar = gsl_matrix_alloc(P, P);
	gsl_multifit_covar(s->J, 0.0, covar);

#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	double chi = gsl_blas_dnrm2(s->f);
	double dof = N - P;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));

	std::cout << "roop(iter): " << iter << std::endl;
	std::cout << "/*********************solution*********************/" << std::endl;
	std::cout << "a	      = " << this->prm_a[1] << " +/- " << c*ERR(0) << std::endl;
	std::cout << "mu1     = " << this->prm_mu1[0] << " +/- " << c*ERR(1) << std::endl;
	std::cout << "mu2     = " << this->prm_mu2[0] << " +/- " << c*ERR(2) << std::endl;
	std::cout << "sigma11 = " << this->prm_sigma1[0] << " +/- " << c*ERR(3) << std::endl;
	std::cout << "sigma11 = " << this->prm_sigma2[0] << " +/- " << c*ERR(4) << std::endl;
	std::cout << "sigma21 = " << this->prm_sigma1[1] << " +/- " << c*ERR(5) << std::endl;
	std::cout << "sigma22 = " << this->prm_sigma2[1] << " +/- " << c*ERR(6) << std::endl;
	std::cout << "k1      = " << this->prm_kappa[0] << " +/- " << c*ERR(7) << std::endl;
	std::cout << "k2      = " << this->prm_kappa[1] << " +/- " << c*ERR(8) << std::endl;
	std::cout << "k3      = " << this->prm_kappa[2] << " +/- " << c*ERR(9)<< std::endl;
	std::cout << "/**************************************************/" << std::endl;
	std::cout << "status  = " << gsl_strerror(status) << "\n" << std::endl;

	gsl_multifit_fdfsolver_free(s);

	return gsl_strerror(status);
}

// 変位のpdfを求める
void Analysis::createDispPdf()
{
	std::cout << "Creating a file of the displacement pdf(.dat).\n" << std::endl;

	// カウント変数
	size_t tmp;

	double gsa_xmin	= -6., integration = 0.;

	FILE *gsax1pdf;
	gsax1pdf = fopen("gsay1pdf.dat", "w");

	for (tmp = 0; tmp*0.01 < 12; tmp++)
	{
		fprintf(gsax1pdf, "%lf %lf\n", gsa_xmin + tmp*0.01, this->_createGaussianPdf(this->prm_a, this->prm_mu1, this->prm_sigma1, gsa_xmin + tmp*0.01));
		integration += 0.01*this->_createGaussianPdf(this->prm_a, this->prm_mu1, this->prm_sigma1, gsa_xmin + tmp*0.01);
	}

	fclose(gsax1pdf);
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
void Analysis::createLevelCrossing()
{
	std::cout << "Creating a file of the level crossing(.dat).\n" << std::endl;

	// カウント変数
	size_t tmp_xi;

	FILE *x1_prob_pass;
	x1_prob_pass = fopen("x1_prob_pass.dat", "w");

	double pp_xi, prob_pass = 0.;

	for (tmp_xi = 0; tmp_xi*0.01 < 8; tmp_xi++)
	{
		prob_pass = 0.;	// 元に戻す
		pp_xi = (double)tmp_xi*0.01;

		prob_pass = this->_culcLevelCrossing(pp_xi, this->prm_a, this->prm_mu1, this->prm_mu2, this->prm_sigma1, this->prm_sigma2, this->prm_kappa);

		fprintf(x1_prob_pass, "%lf %lf\n", pp_xi, prob_pass);
	}

	fclose(x1_prob_pass);
}

// 速度のpdfを作成する
void Analysis::createVelPdf()
{
	std::cout << "Creating a file of the velocity varience(.dat).\n" << std::endl;

	// カウント変数
	size_t tmp;

	double gsa_ymin = -15., integration = 0.;

	FILE *gsax2pdf;
	gsax2pdf = fopen("gsay2pdf.dat", "w");

	for (tmp = 0; tmp*0.01 < 30; tmp++)
	{
		fprintf(gsax2pdf, "%lf %lf\n", gsa_ymin + tmp*0.01, this->_createGaussianPdf(this->prm_a, this->prm_mu2, this->prm_sigma2, gsa_ymin));
		integration += 0.01*this->_createGaussianPdf(this->prm_a, this->prm_mu2, this->prm_sigma2, gsa_ymin);
	}

	fclose(gsax2pdf);
	printf("y integration = %lf \n",integration);
}

/**
 * 等価線形化法により初期値を計算
 *
 */
void Analysis::_culcInitValue(double lambda, double beta2, double ggd_kappa, double alpha, double *sigma_x, double *sigma_y, double *rho_xy)
{
	int tmp, s;
	
	double a[iN*iN], b[iN];					// 線形連立方程式Ax=bのAとB
	double Exxold, Exyold, Eyyold, Exx=0., Exy=0., Eyy=0., err;	// ループ計算時の解の保存用，x:変位，y:速度

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a	= sqrt(gsl_sf_gamma(1./ggd_kappa)*pow(gsl_sf_gamma(3./ggd_kappa), -1.)*beta2);


	/******* 初期値 *******
	 * a_1 = 1, a_2 = a_3 = 0
	 * mu_1 = mu_2 = 0：対称なガウス分布を仮定しているから（ただ，後で mu_1 = 1 にするかも）
	 * sigma_1, sigma_2 はそのまま
	 * kappa_1 = 0（ro をすべて0で仮定し，平均もゼロだから）
	 */

	int loop = 0;
	/******* 応答のガウス性を仮定し，等価線形化法により2次定常モーメントを求める *******/
	do{
		double ke = 1.+3.*EPSILON*Exx;	// 等価線形係数 ke = 1+3εE[X^2]

		for(tmp=0; tmp<iN*iN; tmp++) a[tmp] = 0.;
		for(tmp=0; tmp<iN; tmp++) b[tmp]	= 0.;
		
		/**************** 定常モーメント方程式の係数をAとBに代入 ********************/
		// dXX
		a[0*iN+1] = 2.;	b[0] = 0.;
				
		// dXY
		a[1*iN+0] = -ke;		a[1*iN+1] = -2.*ZETA;	a[1*iN+2] = 1.;	b[1] = 0.;
		
		// dYY
		a[2*iN+1] = 2.*ke;	a[2*iN+2] = 4.*ZETA;	b[2] = alpha*2.*PI*S0 + (1.-alpha)*lambda*(pow(ggd_a,2.)*gsl_sf_gamma(3./ggd_kappa)*pow(gsl_sf_gamma(1./ggd_kappa),-1.));
		/**************************************************************************/
		
		gsl_matrix_view m	= gsl_matrix_view_array(a, iN, iN);
		gsl_vector_view c	= gsl_vector_view_array(b, iN);
		gsl_vector *x		= gsl_vector_alloc(iN);
		gsl_permutation *px	= gsl_permutation_alloc(iN);
		
		// LU分解の方法でモーメント方程式を解く
		gsl_linalg_LU_decomp(&m.matrix, px, &s);
		gsl_linalg_LU_solve(&m.matrix, px, &c.vector, x);
		
		Exx = gsl_vector_get(x,0);
		Exy = gsl_vector_get(x,1);
		Eyy = gsl_vector_get(x,2);
		
		gsl_vector_free(x);
		gsl_permutation_free(px);

		Exxold = Exx;
		Exyold = Exy;
		Eyyold = Eyy;
		loop += 1;
		if (loop == 0) continue;

		err = pow(Exx - Exxold, 2.) + pow(Exy - Exyold, 2.) + pow(Eyy - Eyyold, 2.);	// 収束条件に使う誤差，前ループとの差の二乗和

	} while(err > 10e-6);

	*sigma_x = sqrt(Exx);
	*sigma_y = sqrt(Exy);
	*rho_xy  = Exy/sqrt(Exx*Eyy);

	std::cout << "sigma_x=" << sqrt(Exx) << ", sigma_y=" << sqrt(Eyy) << ", rho_xy=" << Exy/sqrt(Exx*Eyy) << "\n" << std::endl;
	return;
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
double Analysis::_createGaussianPdf(double a[], double mu[], double sigma[], double x) 
{
	double pdf = 0.;
	size_t tmp;

	for (tmp = 0; tmp < NUM_GAUSS; tmp++)
	{
		pdf += a[tmp]*(1./sqrt(2.*PI)/sigma[tmp])*exp(-(x-mu[tmp])*(x-mu[tmp])/2./pow(sigma[tmp],2.));
	}

	return pdf;
}

/**
 * 閾値通過率を計算
 *
 * 
 */
double Analysis::_culcLevelCrossing(double pp_xi, double a[], double mu1[], double mu2[], double sigma1[], double sigma2[], double kappa[])
{
	double prob_pass;
	double pp_c, pp_g, pp_sigma;
	size_t tmp;

	for (tmp = 0; tmp < 3; tmp++)
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

int main(int argc, char *argv[])
{
	Analysis *ana;
	ana = new Analysis;

	// 入力パラメータ
	char *ends;
	ana->lambda = strtod(argv[1], &ends);
	ana->beta2 = strtod(argv[2], &ends);
	ana->alpha = strtod(argv[3], &ends);
	ana->mu1 = strtod(argv[4], &ends);
	ana->mu2 = ana->mu1;

	std::cout << "--------------------------\n" << std::endl;
	std::cout << "lsm.cpp started.\n" << std::endl;

	std::string result = ana->leastSquareMethod();
	if (result == "success")
	{
		ana->createDispPdf();
	}

	std::cout << "lsm.cpp has done.\n" << std::endl;
	std::cout << "--------------------------\n" << std::endl;

	return 0;
}