#ifndef __GA_H_INCLUDE__
#define __GA_H_INCLUDE__

#include <stdlib.h>
#include <random>
#include <math.h>
#define ARRAY_LENGTH(array) (sizeof(array) / sizeof(array[0]))

class GA
{
public:
	GA();
	~GA();
	void initGene();
	void culcFitness();
	void output(int generation);

private:
	/*
	適応度の計算指標
		0:未使用
		1:適応度計算前（突然変異はこの個体だけに適応）
		2:適応度計算済み（交叉時に親とみなす）
	*/
	int* fitnessIndex;
	double* fitness;
	double meanFitness;
	double maxFitness;
	int population;
	int geneLength;
	unsigned int **allIndividual;
	double _binary2Decimal(int* binary);
	bool _isDuplicatedGene(unsigned int **gene);
	// @TODO:setObjectiveFuncを追加したり
	double _getObjectiveFunc(double x);	// 1変数の場合
};

#endif // !__GA_H_INCLUDE__