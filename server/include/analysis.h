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

class GAIndividual;

class Analysis
{
private:
	static const double GGD_KAPPA;
	double _lambda, _beta2, _alpha, _mu1, _mu2;
	void _culcInitValue(double *sigma_x, double *sigma_y, double *rho_xy);

public:
	explicit Analysis(double, double, double, double, double);
	std::string leastSquareMethod(std::vector<double> &);
	int GeneticAlgorithm(std::vector<GAIndividual> &);
	void outputPopsIntoFile(const std::string, const GAIndividual &, const std::vector<double> &, const std::vector<double> &);
	void outputAllPopsIntoFile(const std::string, const std::vector<GAIndividual> &);
	~Analysis();
};

#endif // !__ANALISYS_H_INCLUDE__