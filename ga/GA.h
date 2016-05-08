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
	void culcFitness();
	void outputGeneration(int generation);
	void uniformCrossover();
	int selectIndividual();
	void mutation(double mutationRate);
	void selectRanking();

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

	/*** NSGA-2用 ***/
	std::vector< std::vector<int> > _searchPopulation;
	std::vector< std::vector<int> > _archivePopulation;
	void _initSearchPopulation();
	void _joinGene(const std::vector <std::vector<int> >, const std::vector <std::vector<int> >, std::vector <std::vector<int> >);
	void _nonSuperioritySort(const std::vector <std::vector<int> >, std::vector<std::vector<std::vector<int> > >);
	void _binary2ObjectiveFunc(const std::vector<int>, std::vector<double>);
};

// テンプレート関数の実装を記述したファイルを読み込む
// link:http://qiita.com/MasayaMizuhara/items/37f8c2a5462a4f7f8ea0
#include "detail/GA.h"

#endif // !__GA_H_INCLUDE__