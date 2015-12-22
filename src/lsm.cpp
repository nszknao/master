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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include "expfit.cpp"

#define N 15
#define P 10
#define PI M_PI
#define iN 3	// 初期値を求めるときの方程式の数

/************ 系の係数・入力条件（不変）*******************/
#define S0 (1./(2.*PI))
#define zeta 0.2
#define epi 0.3
#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布


void print_state (size_t iter, gsl_multifit_fdfsolver * s);
void init_values(double lambda, double beta2, double ggd_kappa, double alpha, double *sigma_x, double *sigma_y, double *rho_xy);
double create_gaussian_pdf(double a[3], double mu[3], double sigma[3], double x);
double culc_level_crossing(double zeta, double a[3], double mu1[3], double mu2[3], double sigma1[3], double sigma2[3], double kappa[3]);

int main (int argc, char *argv[])
{
	// カウント変数
	size_t i, tmp;

	char *ends;
	double lambda	= strtod(argv[1], &ends);
	double beta2	= strtod(argv[2], &ends);
	double alpha	= strtod(argv[3], &ends);

	double mu1		= strtod(argv[4], &ends);
	double mu2		= mu1;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a	= sqrt(gsl_sf_gamma(1./GGD_KAPPA)*pow(gsl_sf_gamma(3./GGD_KAPPA), -1.)*beta2);

	// 入力に関するモーメント
	double dF[6];
	dF[0]	= 0;
	dF[1]	= alpha*S0+lambda*(1.-alpha)*(pow(ggd_a,2.)*gsl_sf_gamma(3./GGD_KAPPA)*pow(gsl_sf_gamma(1./GGD_KAPPA),-1.));
	dF[2]	= 0;
	dF[3]	= lambda*pow((1.-alpha),2.)*(pow(ggd_a,4.)*gsl_sf_gamma(5./GGD_KAPPA)*pow(gsl_sf_gamma(1./GGD_KAPPA),-1.));
	dF[4]	= 0;
	dF[5]	= lambda*pow((1.-alpha),3.)*(pow(ggd_a,6.)*gsl_sf_gamma(7./GGD_KAPPA)*pow(gsl_sf_gamma(1./GGD_KAPPA),-1.));

	const size_t n	= N;
	const size_t p	= P;

	/* 近似対象となる観測データを生成 */
	double y[N];
	for (tmp = 0; tmp < N; tmp++) y[tmp] = 0.;

	struct data d = {n, p, y, zeta, epi, dF};
	
	// 乱数生成器
	gsl_rng * r;
	gsl_vector_view x;
	gsl_rng_env_setup();
	const gsl_rng_type * type;
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

	gsl_multifit_function_fdf f;
	f.f			= &expb_f;
	f.df		= &expb_df;
	f.fdf		= &expb_fdf;
	f.n			= n;
	f.p			= p;
	f.params	= &d;
	
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	
	// 等価線形化法で初期値を計算
	double sigma_x, sigma_y, rho_xy;
	init_values (lambda, beta2, GGD_KAPPA, alpha, &sigma_x, &sigma_y, &rho_xy);

	/*  初期値         {a,   μ1,           μ2,          σ11,     σ12,     σ21,     σ22,     k1,                    k2, k3}*/
	double x_init[P] = {0.5, sigma_x+mu1, sigma_y+mu2, sigma_x, sigma_y, sigma_x, sigma_y, rho_xy*sigma_x*sigma_y, 0., 0.};

	//**********************************************************************/
	x = gsl_vector_view_array (x_init, p);
	gsl_multifit_fdfsolver_set(s, &f, &x.vector);
	size_t iter = 0;
	int status;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if (status) break;
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-10, 1e-10);
	} while(status == GSL_CONTINUE && iter < 10000);

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar(s->J, 0.0, covar);
	
	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));

	printf("\nroop(iter): %d\n", iter) ;
	printf("/*********************solution*********************/\n");
	printf("a = %.5lf +/- %.5lf\n",   FIT(0), c*ERR(0));
	printf("mu1 = %.5lf +/- %.5lf\n", FIT(1), c*ERR(1));
	printf("mu2 = %.5lf +/- %.5lf\n", FIT(2), c*ERR(2));
	printf("sigma11 = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
	printf("sigma12 = %.5f +/- %.5f\n", FIT(4), c*ERR(4));
	printf("sigma21 = %.5f +/- %.5f\n", FIT(5), c*ERR(5));
	printf("sigma22 = %.5f +/- %.5f\n", FIT(6), c*ERR(6));
	printf("k1 = %.5f +/- %.5f\n", FIT(7), c*ERR(7));
	printf("k2 = %.5f +/- %.5f\n", FIT(8), c*ERR(8));
	printf("k3 = %.5f +/- %.5f\n", FIT(9), c*ERR(9));
	printf("/**************************************************/\n");
	printf("status = %s\n", gsl_strerror (status));


	/********** パラメータ **********/
	double prm_a[3], prm_mu1[3], prm_mu2[3], prm_sigma1[3], prm_sigma2[3], prm_kappa[3];
	// 重み
	prm_a[1]	= FIT(0);				// a2
	prm_a[0]	= (1.-prm_a[1])/2.;	// a1
	prm_a[2]	= prm_a[0];			// a3
	// 変位
	prm_mu1[0]	= FIT(1);			// mu11
	prm_mu1[1]	= 0.;				// mu21
	prm_mu1[2]	= -1.*prm_mu1[0];	// mu31
	prm_sigma1[0]	= FIT(3);			// sigma11
	prm_sigma1[1]	= FIT(5);			// sigma21
	prm_sigma1[2]	= prm_sigma1[0];	// sigma31
	// 速度
	prm_mu2[0]	= FIT(2);			// mu12
	prm_mu2[1]	= 0.;				// mu22
	prm_mu2[2]	= -1.*prm_mu2[0];	// mu32
	prm_sigma2[0]	= FIT(4);			// sigma12
	prm_sigma2[1]	= FIT(6);			// sigma22
	prm_sigma2[2]	= prm_sigma2[0];	// sigma32
	// 共分散
	prm_kappa[0]	= FIT(7);			// kappa123
	prm_kappa[1]	= FIT(8);			// kappa223
	prm_kappa[2]	= FIT(9);			// kappa323


	if (strcmp(gsl_strerror (status), "success") == 0) {

		/********** 変位応答のpdfを生成 **********/
		#define gsa_xmin -6.
		FILE *gsax1pdf;
		gsax1pdf	= fopen("gsay1pdf.dat","w");
		double integration = 0.;
		for(tmp = 0; tmp*0.01 < 12; tmp++)
		{
			fprintf(gsax1pdf, "%lf %lf\n", gsa_xmin + tmp*0.01, create_gaussian_pdf(prm_a, prm_mu1, prm_sigma1, gsa_xmin + tmp*0.01));
			integration += 0.01*create_gaussian_pdf(prm_a, prm_mu1, prm_sigma1, gsa_xmin + tmp*0.01);
		}
		fclose(gsax1pdf);
//		printf("y1 integration = %lf \n",integration);
		/********** 変位分散を計算 **********/
		FILE *y1_var;
		y1_var  = fopen("anl_y1_var.dat", "w");
		gsax1pdf = fopen("gsay1pdf.dat", "r");
		
		double var_y1 = 0.;
		double mo4_y1 = 0.;
		int line = 0, ret1;
		double row1, row2;
		while((ret1 = fscanf(gsax1pdf, "%lf %lf", &row1, &row2)) != EOF)
		{
			if((line%10) == 0)
			{
				var_y1 += pow(row1,2.)*row2;
				mo4_y1 += pow(row1,4.)*row2;
			}
			line++;
		}
		
		fprintf(y1_var, "%lf", mo4_y1/pow(var_y1,2.));
	
		fclose(gsax1pdf);
		fclose(y1_var);

		/********** 閾値通過率を計算 **********/
		FILE *x1_prob_pass;
		x1_prob_pass	= fopen("x1_prob_pass.dat", "w");

		size_t tmp_xi;
		double pp_xi, prob_pass = 0.;

		for (tmp_xi = 0; tmp_xi*0.01 < 8; tmp_xi++)
		{
			prob_pass	= 0.;	// 元に戻す
			pp_xi	= (double)tmp_xi*0.01;

			prob_pass	= culc_level_crossing(pp_xi, prm_a, prm_mu1, prm_mu2, prm_sigma1, prm_sigma2, prm_kappa);

			fprintf(x1_prob_pass, "%lf %lf\n", pp_xi, prob_pass);
		}

		fclose(x1_prob_pass);

		/********** 速度応答のpdfを生成 **********/
//		#define gsa_ymin -15.
//		FILE *gsax2pdf;
//		gsax2pdf=fopen("gsay2pdf.dat","w");
//		integration = 0.;
//		for(tmp = 0; tmp*0.01 < 30; tmp++) 
//		{
//			fprintf(gsax2pdf, "%lf %lf\n", gsa_ymin+tmp*0.01, create_gaussian_pdf(prm_a, prm_mu2, prm_sigma2, gsa_ymin));
//			integration += 0.01*create_gaussian_pdf(prm_a, prm_mu2, prm_sigma2, gsa_ymin);
//		}
//		fclose(gsax2pdf);
////	printf("y integration = %lf \n",integration);
	}

	gsl_multifit_fdfsolver_free(s);
	return 0;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	printf("iter: %3u x = % 15.8f % 15.8f  |f(x)| = %g\n",
	iter,
	gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	gsl_blas_dnrm2 (s->f));
}

void init_values(double lambda, double beta2, double ggd_kappa, double alpha, double *sigma_x, double *sigma_y, double *rho_xy)
{
	int tmp, s;
	
	double a[iN*iN], b[iN];					// 線形連立方程式Ax=bのAとB
	double Exxold, Exyold, Eyyold, A=0., B=0., C=0., err;	// ループ計算時の解の保存用，x:変位，y:速度

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a		= sqrt(gsl_sf_gamma(1./ggd_kappa)*pow(gsl_sf_gamma(3./ggd_kappa), -1.)*beta2);


	/******* 初期値 *******
	 * a_1 = 1, a_2 = a_3 = 0
	 * mu_1 = mu_2 = 0：対称なガウス分布を仮定しているから（ただ，後で mu_1 = 1 にするかも）
	 * sigma_1, sigma_2 はそのまま
	 * kappa_1 = 0（ro をすべて０で仮定し，平均もゼロだから）
	 */

	/******* 応答のガウス性を仮定し，等価線形化法により2次定常モーメントを求める *******/
	do{
		Exxold = A;
		Exyold = B;
		Eyyold = C;
		
		double ke = 1.+3.*epi*A;	// 等価線形係数 ke = 1+3εE[X^2]

		for(tmp=0; tmp<iN*iN; tmp++) a[tmp] = 0.;
		for(tmp=0; tmp<iN;   tmp++) b[tmp] = 0.;
		
		/**************** 定常モーメント方程式の係数をAとBに代入 ********************/
		// dXX
		a[0*iN+1] = 2.;	b[0] = 0.;
				
		// dXY
		a[1*iN+0] = -ke;		a[1*iN+1] = -2.*zeta;	a[1*iN+2] = 1.;	b[1] = 0.;
		
		// dYY
		a[2*iN+1] = 2.*ke;	a[2*iN+2] = 4.*zeta;	b[2] = alpha*2.*PI*S0 + (1.-alpha)*lambda*(pow(ggd_a,2.)*gsl_sf_gamma(3./ggd_kappa)*pow(gsl_sf_gamma(1./ggd_kappa),-1.));
		/**************************************************************************/
		
		gsl_matrix_view m	= gsl_matrix_view_array(a, iN, iN);
		gsl_vector_view c	= gsl_vector_view_array(b, iN);
		gsl_vector *x		= gsl_vector_alloc(iN);
		gsl_permutation *px	= gsl_permutation_alloc(iN);
		
		// LU分解の方法でモーメント方程式を解く
		gsl_linalg_LU_decomp(&m.matrix, px, &s);
		gsl_linalg_LU_solve(&m.matrix, px, &c.vector, x);
		
		A = gsl_vector_get(x,0);
		B = gsl_vector_get(x,1);
		C = gsl_vector_get(x,2);
		err = pow(A-Exxold,2.)+pow(B-Exyold,2.)+pow(C-Eyyold,2.);	// 収束条件に使う誤差，前ループとの差の二乗和
		
		gsl_vector_free(x);
		gsl_permutation_free(px);
		
	} while(err > 10e-6);

	*sigma_x = sqrt(A);
	*sigma_y = sqrt(C);
	*rho_xy  = B/sqrt(A*C);

	printf("sigma_x=%lf sigma_y=%lf rho_xy=%lf\n", sqrt(A), sqrt(C), B/sqrt(A*C));
	return;
}

// 1変数のガウス分布を3つ足し合わせる関数
double create_gaussian_pdf(double a[3], double mu[3], double sigma[3], double x) 
{
	double gaussian_pdf = 0.;
	size_t tmp;

	for (tmp = 0; tmp < 3; tmp++)
	{
		gaussian_pdf += a[tmp]*(1./sqrt(2.*PI)/sigma[tmp])*exp(-(x-mu[tmp])*(x-mu[tmp])/2./pow(sigma[tmp],2.));
	}

	return gaussian_pdf;
}

// 閾値通過率を計算
double culc_level_crossing(double zeta, double a[3], double mu1[3], double mu2[3], double sigma1[3], double sigma2[3], double kappa[3])
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