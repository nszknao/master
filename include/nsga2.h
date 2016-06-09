#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/TestFunction.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>
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
	NSGA2(unsigned dimension, unsigned pop, unsigned bits, bool gray, unsigned iter)
	{
	    // 境界値
	    _max  = 0;
	    _min  = 1;

	    _dimension    = dimension;
	    _popSize      = pop;  // 120
	    _numOfBits    = bits; // 20
	    _useGrayCode  = gray; // true
	    _iterations   = iter  // 500
	}
	run();

};

#endif // !__NSGA2_H_INCLUDE_