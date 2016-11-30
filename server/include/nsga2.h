#ifndef __NSGA2_H_INCLUDE_
#define __NSGA2_H_INCLUDE_

#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

using namespace std;

class GAIndividual;
class MomentEq;

class NSGA2
{
private:
	std::vector<GAIndividual> _finalPops;
	std::size_t _popSize, _iterations;
	ArchiveMOO _archive;

	void _setValueRange(std::vector<double> &, std::vector<double> &);
	void _saveArchive(ArchiveMOO &);

public:
	explicit NSGA2(std::size_t, std::size_t);
	int run(MomentEq *);
	std::vector<GAIndividual> getFinalPops();
	~NSGA2();
};

#endif // !__NSGA2_H_INCLUDE_
