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

private:
	/*
		適応度の計算指標
			0:未使用
			1:適応度計算前（突然変異はこの個体だけに適応）
			2:適応度計算済み（交叉時に親とみなす）
		@TODO:enumを用意する
	*/
	std::vector<int> fitnessIndex;
	std::vector<double> fitness;
	double meanFitness;
	double maxFitness;
	int maxFitnessNumber;
	int population;
	int geneLength;
	std::vector<std::vector<int>> allIndividual;
	double _binary2Phenotype(int* binary);
	bool _isDuplicatedGene(unsigned int **gene);
	// @TODO:setObjectiveFuncを追加したり
	double _getObjectiveFunc(double x);	// 1変数の場合
};

#endif // !__GA_H_INCLUDE__