#include "../GaCommon.h"
#include <vector>
#include <iostream>

/*
	要素で指定して2次元配列の要素を削除
	削除された要素のインデックスは詰められる
	同じ要素はすべて削除される
*/
template<typename T>
void GaCommonTemp<T>::removeElement(std::vector<std::vector<T> > &population, const std::vector<T> &element)
{
	for (int tmp = 0; tmp < population.size(); ++tmp)
		if (population[tmp] == element) population.erase(population.begin() + tmp);
}

/*
	2次元配列のすべての要素を表示する
	@param &individual 表示する2次元配列
*/
template<typename T>
void GaCommonTemp<T>::outputAllElement(const std::vector<T> &individual)
{
	for (int tmp = 0; tmp < individual.size(); ++tmp)
	{
		std::cout << individual[tmp];
		if (tmp != individual.size()-1)
			std::cout << " ";
	}
	std::cout << std::endl;
}