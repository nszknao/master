#ifndef __GAINDIVIDUAL_H_INCLUDE_
#define __GAINDIVIDUAL_H_INCLUDE_

#include <vector>

struct GAIndividual
{
	int index;
	std::vector<double> pValue;
	std::vector<double> oValue;
	std::vector<double> mValue;
};

#endif // !__GAINDIVIDUAL_H_INCLUDE_