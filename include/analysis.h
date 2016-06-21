#ifndef __ANALISYS_H_INCLUDE__
#define __ANALISYS_H_INCLUDE__

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

#define NUM_OF_MOMENTEQ 15
#define NUM_OF_PARAM 10
#define PI M_PI
#define NUM_GAUSS 3	// 足しあわせるガウス分布の数

/************ 系の係数・入力条件（不変）*******************/
#define S0 1./(2.*PI)
#define ZETA 0.2
#define EPSILON 0.3
#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布


class Analysis
{
private:
	double _lambda, _beta2, _alpha, mu1, mu2;
	// 混合ガウス分布のパラメータ
	double prm_a[NUM_GAUSS], prm_mu1[NUM_GAUSS], prm_mu2[NUM_GAUSS], prm_sigma1[NUM_GAUSS], prm_sigma2[NUM_GAUSS], prm_kappa[NUM_GAUSS];
	double _createGaussianPdf(double a[], double mu[], double sigma[], double x);
	double _culcLevelCrossing(double pp_xi, double a[], double mu1[], double mu2[], double sigma1[], double sigma2[], double kappa[]);
	void _culcInitValue(double *sigma_x, double *sigma_y, double *rho_xy);

public:
	// TODO:グローバル変数はアンスコから始める．
	Analysis(double lambda, double beta2, double alpha, double mu1, double mu2);
	std::string leastSquareMethod();
	void createDispPdf();
	void createVelPdf();
	void culcDispVarience();
	void culcVelVarience();
	void createLevelCrossing();
};

#endif // !__ANALISYS_H_INCLUDE_