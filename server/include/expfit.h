#ifndef __EXPFIT_H_INCLUDE_
#define __EXPFIT_H_INCLUDE_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

class MomentEq
{
public:
	static std::vector<double> expb_f(const std::vector<double> &, const std::vector<double> &dG, const std::vector<std::size_t> &f);
	static void getMomentFromParameter(const std::vector<double> &, std::vector<double> &);
};

#endif // !__EXPFIT_H_INCLUDE_
