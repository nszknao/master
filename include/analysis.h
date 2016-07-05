#ifndef __ANALISYS_H_INCLUDE__
#define __ANALISYS_H_INCLUDE__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
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
#define NUM_GAUSS 3	// 足しあわせるガウス分布の数
#define PI M_PI

/************ 系の係数・入力条件（不変）*******************/
#define S0 1./(2.*PI)
#define ZETA 0.2
#define EPSILON 0.3
#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布


class Parameter
{
private:
	std::vector<double> _a, _mu1, _mu2, _sigma1, _sigma2, _kappa;
public:
	void setParameter(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void allocParameter();
	void freeParameter();
	std::vector<double> getParameter(std::string);
};

class Analysis
{
private:
	double _lambda, _beta2, _alpha, mu1, mu2;
	double _createGaussianPdf(const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, double x);
	double _culcLevelCrossing(double pp_xi, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &);
	void _culcInitValue(double *sigma_x, double *sigma_y, double *rho_xy);

public:
	// TODO:グローバル変数はアンスコから始める．
	Analysis(double lambda, double beta2, double alpha, double mu1, double mu2);
	std::string leastSquareMethod(Parameter *prm);
	void createDispPdf(Parameter *prm);
	void createVelPdf(Parameter *prm);
	void culcDispVarience();
	void culcVelVarience();
	void createLevelCrossing(Parameter *prm);
	~Analysis();
};

#endif // !__ANALISYS_H_INCLUDE_