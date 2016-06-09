#ifndef __MOMENTEQ_H_INCLUDE_
#define __MOMENTEQ_H_INCLUDE_

#include <iostream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "../include/paramdata.h"

using namespace std;

class MomentEqFunc
{
private:
	double *_r;	// 計算結果を格納
public:
	MomentEqFunc(const std::vector<double>& x, ParamData *params);
	double momentEqFunc::F1();
	double momentEqFunc::F2();
	double momentEqFunc::F3();
	double momentEqFunc::F4();
	double momentEqFunc::F5();
	double momentEqFunc::F6();
	double momentEqFunc::F7();
	double momentEqFunc::F8();
	double momentEqFunc::F9();
	double momentEqFunc::F10();
	double momentEqFunc::F11();
	double momentEqFunc::F12();
	double momentEqFunc::F13();
	double momentEqFunc::F14();
	double momentEqFunc::F15();
};

#endif // !__MOMENTEQ_H_INCLUDE_