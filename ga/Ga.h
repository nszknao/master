#ifndef __GA_H_INCLUDE__
#define __GA_H_INCLUDE__

#define _USE_MATH_DEFINES

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
	void initGene();
	void outputGeneration(int generation);
	void uniformCrossover();
	int selectIndividual();
	void mutation(double mutationRate);
	void selectRanking();

	/*** NSGA-2—p ***/
	void nsga2Run();
private:
	std::vector<std::vector<int> > allIndividual;
	std::vector<double> fitness;
	int _numVariable;
	int _population;
	int _geneLength;

	double _binary2Phenotype(const std::vector<int>);
	void _convertPhenotype(const std::vector<int>, std::vector<double>);
	bool _isDuplicatedGene(std::vector <std::vector<int> > &gene, int column);
	void _outputIndividuals(std::vector <std::vector<int> > &individuals);
	double _getObjectiveFunc(const std::vector<double>, std::vector<double>);
	double _f1(const std::vector<double>);
	double _f2(const std::vector<double>);

	/*** NSGA-2—p ***/
	std::vector< std::vector<int> > _searchPopulation;
	std::vector< std::vector<int> > _archivePopulation;
	void _initSearchPopulation();
	void _nonSuperioritySort(const std::vector <std::vector<int> >, std::vector<std::vector<std::vector<int> > >);
	void _binary2ObjectiveFunc(const std::vector<int>, std::vector<double>);
	void _sortByObjectiveValue(const std::vector<std::vector<int> >, std::vector<std::vector<int> >, int num);
	void _updateArchivePopulation(const std::vector<std::vector<std::vector<int> > >, std::vector<std::vector<int> >);
	void _crowdingSort(const std::vector<std::vector<int > >, std::vector<std::vector<int> >);
	double _culcCrowdingDistanse(const std::vector<std::vector<std::vector<int> > >, int, int);
	double _culcCrowdingDistanseForIndividual(const std::vector<std::vector<std::vector<int> > > , const std::vector<int>);
	void _insertIndividuals(std::vector<std::vector<int> >, const std::vector<std::vector<int> >);
	void _crowdedTournamentSlection(const std::vector<std::vector<int> >, std::vector<std::vector<int> >, const std::vector<std::vector<std::vector<int> > >);
	int _returnGeneRank(const std::vector<std::vector<std::vector<int> > >, const std::vector<std::vector<int> >);
};

#endif // !__GA_H_INCLUDE__