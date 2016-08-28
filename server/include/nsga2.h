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
#include "common.h"

using namespace std;

struct GAIndividual
{
	int index;
	std::vector<double> pValue;
	std::vector<double> oValue;
	std::vector<double> mValue;
};

class NSGA2
{
private:
	std::vector<GAIndividual> _finalPops;
	unsigned _dimension, _popSize, _iterations;
	ArchiveMOO _archive;

	void _setValueRange(std::vector<double> &, std::vector<double> &);
	void _saveArchive(ArchiveMOO &);
	void _allocFinalPops(int);

public:
	NSGA2(int pop, int iter);
	~NSGA2();
	int run(ParamData* f);
	void saveArchiveInFile(const std::string);
	std::vector<GAIndividual> getFinalPops();
};

#endif // !__NSGA2_H_INCLUDE_