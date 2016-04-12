#ifndef __GA_H_INCLUDE__
#define __GA_H_INCLUDE__

#include <stdlib.h>
#include <random>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#define ARRAY_LENGTH(array) (sizeof(array) / sizeof(array[0]))

class GA
{
public:
	GA(int population, int geneLength);
	~GA();
	void initGene();
	void culcFitness();
	void outputGeneration(int generation);
	void uniformCrossover();
	int selectIndividual();
	void mutation(double mutationRate);
	void selectRanking();

private:
	std::vector<std::vector<int>> allIndividual;
	std::vector<double> fitness;
	double meanFitness;	// g‚Á‚Ä‚È‚¢
	double maxFitness;	// g‚Á‚Ä‚È‚¢
	int maxFitnessNumber;	// g‚Á‚Ä‚È‚¢
	int _population;
	int _geneLength;

	void _setPopulation(int population);
	void _setGeneLength(int length);
	int _getPopulation();
	int _getGeneLength();
	// @TODO:setObjectiveFunc‚ğ’Ç‰Á‚µ‚½‚è
	double _getObjectiveFunc(double x);	// 1•Ï”‚Ìê‡
	double _binary2Phenotype(std::vector<int>);
	bool _isDuplicatedGene(std::vector<std::vector<int>> gene, int column);
	void _outputIndividuals(std::vector<std::vector<int>> individuals);
	void _outputIndividual(std::vector<int> individual);
};

#endif // !__GA_H_INCLUDE__