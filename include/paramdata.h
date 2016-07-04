#ifndef __PARAMDATA_H_INCLUDE_
#define __PARAMDATA_H_INCLUDE_

using namespace std;

class ParamData
{
public:
	unsigned int n, p;
	double* dG;
	double zeta, epsilon;

	ParamData(unsigned int arg_n, unsigned int arg_p, double arg_zeta, double arg_epsilon, double* arg_dG)
	{
		n		= arg_n;
		p		= arg_p;
		zeta	= arg_zeta;
		epsilon	= arg_epsilon;
		dG	= arg_dG;
	}
};


#endif // !__PARAMDATA_H_INCLUDE_