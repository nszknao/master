#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

using namespace std;

class GAIndividual;

class NSGA2
{
private:
	std::vector<GAIndividual> _finalPops;
	unsigned int  _dimension, _popSize, _iterations;
	ArchiveMOO _archive;

	void _setValueRange(std::vector<double> &, std::vector<double> &);
	void _saveArchive(ArchiveMOO &);
	static gsl_matrix *_correlationMatrix(const gsl_matrix *);
	std::vector<unsigned int> _dimensionReduction(gsl_vector *, gsl_matrix *);

public:
	NSGA2(int pop, int iter);
	~NSGA2();
	int run(double *);
	std::vector<GAIndividual> getFinalPops();
};

#endif // !__NSGA2_H_INCLUDE_