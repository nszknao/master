#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

#include <gsl/gsl_vector.h>

#include "expfit.h"

using namespace std;


class nsga2_function
{
public:
	size_t n;       /* number of functions */
	size_t p;       /* number of independent variables */
	void * params;  /* user parameters */
};

class NSGA2
{
private:
	unsigned _dimension, _popSize, _numOfBits, _iterations;
	bool _useGrayCode;
	int _max, _min;

	void _saveArchiveInFile(char *filename, ArchiveMOO &archive);

public:
	NSGA2(int pop, int bits, bool gray, int iter);
	int run(nsga2_function * f);
};

#endif // !__NSGA2_H_INCLUDE_