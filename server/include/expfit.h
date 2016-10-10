#ifndef __EXPFIT_H_INCLUDE_
#define __EXPFIT_H_INCLUDE_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/********** 系の係数・入力条件**********/
#ifndef __OSCILLATOR_AND_EXCITATION_PARAMETER__
#define __OSCILLATOR_AND_EXCITATION_PARAMETER__
#define PI M_PI
#define S0 1./(2.*PI)
#define EPSILON 0.3
#define ZETA 0.05
#endif // !__OSCILLATOR_AND_EXCITATION_PARAMETER__

/********** 解析条件 **********/
#ifndef __ANALISYS_PARAMETER__
#define __ANALISYS_PARAMETER__
#define NUM_OF_MOMENTEQ 15
#define NUM_OF_MOMENT 21
#define NUM_OF_PARAM 10
#endif // !__ANALISYS_PARAMETER__

class MomentEq
{
public:
	static int expb_f(const gsl_vector *x, void *params, gsl_vector *f);
	static int expb_df(const gsl_vector * x, void *params, gsl_matrix *J);
	static int expb_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);
	static void getMomentFromParameter(const std::vector<double> &, std::vector<double> &);
	static void getJacobyFromParameter(const std::vector<double> &, std::vector< std::vector<double> > &);
};

#endif // !__EXPFIT_H_INCLUDE_