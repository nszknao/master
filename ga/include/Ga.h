#ifndef __GA_H_INCLUDE__
#define __GA_H_INCLUDE__

//#define _USE_MATH_DEFINES
# define M_PI           3.14159265358979323846

#include <stdlib.h>
#include <random>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>

class GA
{
public:
	GA(int numVariable);
	~GA();

	/*** NSGA-2—p ***/
	void nsga2Run();
private:
	int _numVariable;
	int _population;
	int _geneLength;

	double _binary2Phenotype(const std::vector<int>&);
	void _convertPhenotype(const std::vector<int>&, std::vector<double>&);
	bool _isDuplicatedGene(std::vector <std::vector<int> >&, const std::vector<int> &);
	void _outputGene(const std::vector<int>&);
	void _outputPopulation(const std::vector <std::vector<int> >&);
	void _getObjectiveFunc(const std::vector<double>&, std::vector<double>&);
	double _f1(const std::vector<double>&);
	double _f2(const std::vector<double>&);

	/*** NSGA-2—p ***/
	void _initSearchPopulation(std::vector<std::vector<int> >&);
	void _nonSuperioritySort(const std::vector <std::vector<int> >&, std::vector<std::vector<std::vector<int> > >&);
	void _binary2ObjectiveFunc(const std::vector<int>&, std::vector<double>&);
	void _sortByObjectiveValue(const std::vector<std::vector<int> >&, std::vector<std::vector<int> >&, int num);
	void _updateArchivePopulation(const std::vector<std::vector<std::vector<int> > >&, std::vector<std::vector<int> >&, std::vector<std::vector<int> >&);
	void _crowdingSort(const std::vector<std::vector<int > >&, std::vector<std::vector<int> >&);
	void _putObjectiveSortedGeneEveryObjectiveFunc(const std::vector<std::vector<int> >&, std::vector<std::vector<std::vector<int> > >&);
	double _culcCrowdingDistanse(const std::vector<std::vector<std::vector<int> > >&, int, int);
	double _culcCrowdingDistanseForIndividual(const std::vector<std::vector<std::vector<int> > >&, const std::vector<int>&);
	void _insertIndividuals(std::vector<std::vector<int> >&, const std::vector<std::vector<int> >&);
	void _crowdedTournamentSelection(const std::vector<std::vector<int> >&, std::vector<std::vector<int> >&, const std::vector<std::vector<std::vector<int> > >&);
	int _returnGeneRank(const std::vector<std::vector<std::vector<int> > >&, const std::vector<int>&);
	void _select2GenesFromPopulation(const std::vector<std::vector<int> >&, std::vector<int>&, std::vector<int>&);
	void _mutationGene(std::vector<std::vector<int> >&, double);
	void _uniformCrossover(const std::vector<int>&, const std::vector<int>&, std::vector<int>&, std::vector<int>&);
	void _highRankGeneSelection(const std::vector<std::vector<std::vector<int> > >&, const std::vector<std::vector<int> >&, std::vector<std::vector<int> >&, int num);
	void _outputObjectiveValue(std::vector<std::vector<int> >, int generation);
	bool _isSuperior(const std::vector<int>&, const std::vector<std::vector<int> >&);
	void _createRandomlyIndividual(std::vector<int>&);
	int _numOfSuperior(const std::vector<double>&, const std::vector<double>&);
};

#endif // !__GA_H_INCLUDE__