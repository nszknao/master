#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/TestFunction.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

#include <gsl/gsl_vector.h>

#include "../include/momenteq.h"

using namespace std;

class NSGA2
{
private:
	unsigned _dimension, _popSize, _numOfBits, _iterations;
	bool _useGrayCode;
	int _max, _min;

	void _saveArchiveInFile(char *filename, ArchiveMOO &archive);

public:
	// コンストラクタ
	NSGA2(unsigned pop, unsigned bits, bool gray, unsigned iter) {
	    // 境界値
	    _max  = 0;
	    _min  = 1;

	    _popSize      = pop;  // 120
	    _numOfBits    = bits; // 20
	    _useGrayCode  = gray; // true
	    _iterations   = iter  // 500
	}
	run(nsga2_function * f);
};

class nsga2_function
{
public:
	size_t n;       /* number of functions */
	size_t p;       /* number of independent variables */
	void * params;  /* user parameters */
};

#endif // !__NSGA2_H_INCLUDE_