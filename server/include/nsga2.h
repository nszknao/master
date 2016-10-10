#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

#include <gsl/gsl_vector.h>

/********** 解析条件 **********/
#ifndef __ANALISYS_PARAMETER__
#define __ANALISYS_PARAMETER__
#define NUM_OF_MOMENTEQ 15
#define NUM_OF_MOMENT 21
#define NUM_OF_PARAM 10
#endif // !__ANALISYS_PARAMETER__

using namespace std;

class GAIndividual;

class NSGA2
{
private:
	std::vector<GAIndividual> _finalPops;
	unsigned _dimension, _popSize, _iterations;
	ArchiveMOO _archive;

	void _setValueRange(std::vector<double> &, std::vector<double> &);
	void _saveArchive(ArchiveMOO &);

public:
	NSGA2(int pop, int iter);
	~NSGA2();
	int run(double *);
	std::vector<GAIndividual> getFinalPops();
};

#endif // !__NSGA2_H_INCLUDE_