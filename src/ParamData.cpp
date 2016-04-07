#include "ParamData.h"

ParamData::ParamData(int arg_n, int arg_p, double* arg_y, double arg_zeta, double arg_epsilon, double* arg_dG)
{
	this->n = arg_n;
	this->p = arg_p;
	this->y = arg_y;
	this->zeta = arg_zeta;
	this->epsilon = arg_epsilon;
	this->dG = arg_dG;
}