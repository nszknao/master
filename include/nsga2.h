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
#include "paramdata.h"

using namespace std;

class NSGA2
{
private:
	std::vector< std::vector<double> > _obj, _prm;		// 結果を保存する
	unsigned _dimension, _popSize, _numOfBits, _iterations;
	bool _useGrayCode;
	int _max, _min;

	void _setValueRange(std::vector<double> &, std::vector<double> &);
	void _saveArchiveInFile(char*, ArchiveMOO &);
	void _saveArchive(ArchiveMOO &);

public:
	NSGA2(int pop, int bits, bool gray, int iter);
	void freeVector();
	int run(ParamData* f);
	std::vector< std::vector<double> > getObjValue();
	std::vector< std::vector<double> > getPrmValue();
};

#endif // !__NSGA2_H_INCLUDE_