#ifndef __GA_H_INCLUDE__
#define __GA_H_INCLUDE__

#include <stdlib.h>
#include <random>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#define ARRAY_LENGTH(array) (sizeof(array) / sizeof(array[0]))

class GA
{
public:
	GA();
	~GA();
	void initGene();
	void culcFitness();
	void output(int generation);
	void uniformCrossover();
	int selectIndividual();
	void mutation();
	void selectRanking();

private:
	std::vector<std::vector<int>> allIndividual;
	std::vector<double> fitness;
	double meanFitness;	// g‚Á‚Ä‚È‚¢
	double maxFitness;	// g‚Á‚Ä‚È‚¢
	int maxFitnessNumber;	// g‚Á‚Ä‚È‚¢
	int population;
	int geneLength;

	// @TODO:setObjectiveFunc‚ğ’Ç‰Á‚µ‚½‚è
	double _getObjectiveFunc(double x);	// 1•Ï”‚Ìê‡
	double _binary2Phenotype(int* binary);
	bool _isDuplicatedGene(unsigned int **gene);
};

#endif // !__GA_H_INCLUDE__