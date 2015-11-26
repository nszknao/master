#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "expfit.c"

#define N 15
#define P 10
#define PI M_PI
#define iN 3	// 初期値を求めるときの方程式の数

/************ 系の係数・入力条件（不変）*******************/
#define S0 (1./(2.*PI))
#define zeta 0.05
#define epi 0.3


void print_state (size_t iter, gsl_multifit_fdfsolver * s);
void init_values(double lambda, double beta2, double alpha, double *sigma_x, double *sigma_y, double *rho_xy);

int main (int argc, char *argv[])
{
	// カウント変数
	size_t i;

	char *ends;
	double lambda = strtod(argv[1], &ends);
	double beta2  = strtod(argv[2], &ends);
	double alpha  = strtod(argv[3], &ends);

	double mu1_width = strtod(argv[4], &ends);
//	double mu1_width = 0.489;

	double dG[6];
	dG[0] = 0;
	dG[1] = alpha*S0+lambda*pow(sqrt(1-alpha),2.)*beta2;
	dG[2] = 0;
	dG[3] = 3*lambda*pow(sqrt(1-alpha),4.)*pow(beta2,2.);
	dG[4] = 0;
	dG[5] = 15*lambda*pow(sqrt(1-alpha),6.)*pow(beta2,3.);
	
	const size_t n = N;
	const size_t p = P;
	struct data d = {n, y, zeta, epi, dG};
	
	// 乱数生成器
	gsl_rng * r;
	gsl_vector_view x;
	gsl_rng_env_setup();
	const gsl_rng_type * type;
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

	gsl_multifit_function_fdf f;
	f.f	 = &expb_f;
	f.df	 = &expb_df;
	f.fdf	 = &expb_fdf;
	f.n	 = n;
	f.p	 = p;
	f.params = &d;
	
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	
	/* 近似対象となる観測データを生成 */
	double y[N];
	for (i=0;i<N;i++) y[i] = 0.;

	double sigma_x, sigma_y, rho_xy;
	init_values (lambda, beta2, alpha, &sigma_x, &sigma_y, &rho_xy);

	/*  初期値         {a,   μ1,               μ2,               σ11,    σ12,    σ21,    σ22,    k1,                     k2, k3}*/
	double x_init[P] = {0.5, sigma_x+mu1_width, sigma_y+mu1_width, sigma_x, sigma_y, sigma_x, sigma_y, rho_xy*sigma_x*sigma_y, 0., 0.};

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
	
	// 重み
	double a2 =  gsl_vector_get(s->x,0) ;
	double a1 = (1.0-a2)/2.0;
	double a3 = a1;

	// 変位に関するパラメータ
	double mu11 = gsl_vector_get(s->x,1);
	double mu21 = 0.;
	double mu31 = -1.*mu11;
	double sigma11 = gsl_vector_get(s->x,3);
	double sigma21 = gsl_vector_get(s->x,5);
	double sigma31 = sigma11;

	// 速度に関するパラメータ
	double mu12 = gsl_vector_get(s->x,2);
	double mu22 = 0.;
	double mu32 = -1.*mu12;
	double sigma12 = gsl_vector_get(s->x,4);
	double sigma22 = gsl_vector_get(s->x,6);
	double sigma32 = sigma12;
	
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

	// 重み
	a2 = FIT(0);
	a1 = (1.-a2)/2.;
	a3 = a1;

	// 変位に関するパラメータ
	mu11 = FIT(1);
	mu21 = 0.;
	mu31 = -1.*mu11;
	sigma11 = FIT(3);
	sigma21 = FIT(5);
	sigma31 = sigma11;

	// 速度に関するパラメータ
	mu12 = FIT(2);
	mu22 = 0.;
	mu32 = -1.*mu12;
	sigma12 = FIT(4);
	sigma22 = FIT(6);
	sigma32 = sigma12;


	#define gsa_xmin -6
	#define gsa_ymin -15
	if (strcmp(gsl_strerror (status), "success") == 0) {
/********** 変位応答のpdfを生成 **********/
		FILE *gsax1pdf;
		gsax1pdf=fopen("gsay1pdf.dat","w");
		double integration=0.;
		for(i=0;i<1200;i++) 
		{
			fprintf(gsax1pdf, "%lf %lf\n", gsa_xmin+i*0.01,   a1*(1/sqrt(2*PI)/sigma11)*exp(-(gsa_xmin+i*0.01-mu11)*(gsa_xmin+i*0.01-mu11)/2/sigma11/sigma11) 
									+ a2*(1/sqrt(2*PI)/sigma21)*exp(-(gsa_xmin+i*0.01-mu21)*(gsa_xmin+i*0.01-mu21)/2/sigma21/sigma21) 
									+ a3*(1/sqrt(2*PI)/sigma31)*exp(-(gsa_xmin+i*0.01-mu31)*(gsa_xmin+i*0.01-mu31)/2/sigma31/sigma31));
			integration += 0.01*(a1*(1/sqrt(2*PI)/sigma11)*exp(-(gsa_xmin+i*0.01-mu11)*(gsa_xmin+i*0.01-mu11)/2/sigma11/sigma11)
					+ a2*(1/sqrt(2*PI)/sigma21)*exp(-(gsa_xmin+i*0.01-mu21)*(gsa_xmin+i*0.01-mu21)/2/sigma21/sigma21)
					+ a3*(1/sqrt(2*PI)/sigma31)*exp(-(gsa_xmin+i*0.01-mu31)*(gsa_xmin+i*0.01-mu31)/2/sigma31/sigma31));
		}
		fclose(gsax1pdf);
//		printf("y1 integration = %lf \n",integration);
/********** 変位分散を計算 **********/
		FILE *y1_var;
		y1_var  = fopen("anl_y1_var.dat", "w");
		gsax1pdf = fopen("gsay1pdf.dat", "r");
		
		double var_y1 = 0.;
		double mo4_y1 = 0.;
		int line = 0;
		int ret1;
		double row1, row2;
		while((ret1 = fscanf(gsax1pdf, "%lf %lf", &row1, &row2)) != EOF)
		{
			if((line%10) == 0)
			{
				var_y1 += pow(row1,2)*row2;
				mo4_y1 += pow(row1,4)*row2;
			}
			line++;
		}
		
		fprintf(y1_var, "%lf", mo4_y1/pow(var_y1,2));
	
		fclose(gsax1pdf);
		fclose(y1_var);
/********** 速度応答のpdfを生成 **********/
		FILE *gsax2pdf;
		gsax2pdf=fopen("gsay2pdf.dat","w");
		integration = 0.;
		for(i=0;i<3000;i++) 
		{
			fprintf(gsax2pdf, "%lf %lf\n", gsa_ymin+i*0.01,   a1*(1/sqrt(2*PI)/sigma12)*exp(-(gsa_ymin+i*0.01-mu12)*(gsa_ymin+i*0.01-mu12)/2/sigma12/sigma12)
									+ a2*(1/sqrt(2*PI)/sigma22)*exp(-(gsa_ymin+i*0.01-mu22)*(gsa_ymin+i*0.01-mu22)/2/sigma22/sigma22)
									+ a3*(1/sqrt(2*PI)/sigma32)*exp(-(gsa_ymin+i*0.01-mu32)*(gsa_ymin+i*0.01-mu32)/2/sigma32/sigma32));
			integration += 0.01*(a1*(1/sqrt(2*PI)/sigma12)*exp(-(gsa_ymin+i*0.01-mu12)*(gsa_ymin+i*0.01-mu12)/2/sigma12/sigma12)
					+ a2*(1/sqrt(2*PI)/sigma22)*exp(-(gsa_ymin+i*0.01-mu22)*(gsa_ymin+i*0.01-mu22)/2/sigma22/sigma22)
					+ a3*(1/sqrt(2*PI)/sigma32)*exp(-(gsa_ymin+i*0.01-mu32)*(gsa_ymin+i*0.01-mu32)/2/sigma32/sigma32));
		}
//		printf("y integration = %lf \n",integration);
		fclose(gsax2pdf);
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

/******************** 初期値を求める関数 ********************/
void init_values(double lambda, double beta2, double alpha, double *sigma_x, double *sigma_y, double *rho_xy)
{
	int k, s;
	
	double a[iN*iN], b[iN];					// 線形連立方程式Ax=bのAとB
	double Exxold, Exyold, Eyyold, A=0., B=0., C=0., err;	// ループ計算時の解の保存用，x:変位，y:速度
	double ke;						// 等価線形係数

	/******* 初期値 *******
	 * a_1 = 1, a_2 = a_3 = 0
	 * mu_1 = mu_2 = 0：対称なガウス分布を仮定しているから（ただ，後で mu_1 = 1 にするかも）
	 * sigma_1, sigma_2 はそのまま
	 * kappa_1 = 0（ro をすべて０で仮定し，平均もゼロだから）
	 */

	/******* 応答のガウス性を仮定し，等価線形化法により2次定常モーメントを求める *******/
	do{
		Exxold = A;		 // 前ループの解を保存
		Exyold = B;
		Eyyold = C;
		
		/*************** 等価線形係数 **************************************/
		ke = 1.+3.*epi*A;	// 等価線形係数 ke = 1+3εE[X^2]
		
		for(k=0; k<iN*iN; k++) a[k] = 0.;	// 初期化
		for(k=0; k<iN;   k++) b[k] = 0.;	// 初期化
		
		/**************** 定常モーメント方程式の係数をAとBに代入 ********************/
		// dXX
		a[0*iN+1] = 2.;	b[0] = 0.;
				
		// dXY
		a[1*iN+0] = -ke;		a[1*iN+1] = -2.*zeta;	a[1*iN+2] = 1.;	b[1] = 0.;
		
		// dYY
		a[2*iN+1] = 2.*ke;	a[2*iN+2] = 4.*zeta;	b[2] = alpha*2.*PI*S0 + (1-alpha)*lambda*beta2;
		/**************************************************************************/
		
		gsl_matrix_view m	= gsl_matrix_view_array(a, iN, iN);
		gsl_vector_view c	= gsl_vector_view_array(b, iN);
		gsl_vector *x		= gsl_vector_alloc(iN);
		gsl_permutation *px	= gsl_permutation_alloc(iN);
		
		// LU分解の方法でモーメント方程式を解く
		gsl_linalg_LU_decomp(&m.matrix, px, &s);
		gsl_linalg_LU_solve(&m.matrix, px, &c.vector, x);
		
		A = gsl_vector_get(x,0); // 解の保存
		B = gsl_vector_get(x,1);
		C = gsl_vector_get(x,2);
		err = pow(A-Exxold,2.)+pow(B-Exyold,2.)+pow(C-Eyyold,2.);	// 収束条件に使う誤差，前ループとの差の二乗和
		
		gsl_vector_free(x);
		gsl_permutation_free(px);
		
	}while(err>10e-6);			// 収束判定

	*sigma_x = sqrt(A);
	*sigma_y = sqrt(C);
	*rho_xy  = B/sqrt(A*C);

	return printf("sigma_x=%lf sigma_y=%lf rho_xy=%lf\n", sqrt(A), sqrt(C), B/sqrt(A*C));
}