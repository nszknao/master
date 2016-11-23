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
#include <gsl/gsl_sf.h>

class GAIndividual;

class Analysis
{
private:
	static const double GGD_KAPPA;
	double _lambda, _beta2, _alpha;

public:
	explicit Analysis(double, double, double);
	int GeneticAlgorithm(std::vector<GAIndividual> &);
	void outputAllPopsIntoFile(const std::string, const std::vector<GAIndividual> &);
	~Analysis();
};

#endif // !__ANALISYS_H_INCLUDE__
