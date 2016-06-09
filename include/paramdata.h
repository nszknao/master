#ifndef __PARAMDATA_H_INCLUDE_
#define __PARAMDATA_H_INCLUDE_

using namespace std;

class ParamData
{
public:
	int n, p;
	double* y;
	double* dG;
	double zeta, epsilon;

	ParamData(int arg_n, int arg_p, double* arg_y, double arg_zeta, double arg_epsilon, double* arg_dG)
	{
		n		= arg_n;
		p		= arg_p;
		y		= arg_y;
		zeta	= arg_zeta;
		epsilon	= arg_epsilon;
		dG	= arg_dG;
	}
};


#endif // !__PARAMDATA_H_INCLUDE_