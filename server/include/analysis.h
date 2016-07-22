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

#ifndef __OSCILLATOR_AND_EXCITATION_PARAMETER__
#define __OSCILLATOR_AND_EXCITATION_PARAMETER__
/********** 系の係数・入力条件（不変）**********/
#define PI M_PI
#define S0 1./(2.*PI)
#define EPSILON 0.3
#define ZETA 0.2
#endif // !__OSCILLATOR_AND_EXCITATION_PARAMETER__

/********** 系の係数・入力条件（不変）**********/
#define NUM_OF_MOMENTEQ 15
#define NUM_OF_PARAM 10
#define NUM_GAUSS 3	// 足しあわせるガウス分布の数
#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布


class Parameter
{
private:
	std::vector<double> _a, _mu1, _mu2, _sigma1, _sigma2, _kappa;
public:
	void setParameter(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	void allocParameter();
	void freeParameter();
	bool validate();
	std::vector<double> getParameter(std::string);
};

class Analysis
{
private:
	double _lambda, _beta2, _alpha, _mu1, _mu2;
	double _createGaussianPdf(const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, double x);
	double _culcLevelCrossing(double, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &);
	void _culcInitValue(double *sigma_x, double *sigma_y, double *rho_xy);

public:
	Analysis(double, double, double, double, double);
	std::string leastSquareMethod(Parameter*);
	int GeneticAlgorithm(std::vector<Parameter*> &, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &);
	void createDispPdf(Parameter*, std::vector<double> &, std::vector<double> &, int);
	void createVelPdf(Parameter*, std::vector<double> &, std::vector<double> &, int);
	void createLevelCrossing(Parameter*, std::vector<double> &, std::vector<double> &, int);
	void outputIntoFile(const std::string, const std::vector<double> &, const std::vector<double> &);
	bool isOverSpecifyValue(const std::vector<double> &, double);
	~Analysis();
};

#endif // !__ANALISYS_H_INCLUDE__