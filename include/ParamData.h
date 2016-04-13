#ifndef __PARAMDATA_H_INCLUDE__
#define __PARAMDATA_H_INCLUDE__

class ParamData
{
public:
	int n;
	int p;
	double* y;
	double zeta;
	double epsilon;
	double* dG;

	ParamData(int arg_n, int arg_p, double* arg_y, double arg_zeta, double arg_epsilon, double* arg_dG);
};

#endif