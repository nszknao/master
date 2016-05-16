#include "../include/GaCommon.h"

/*
	２つの二次元配列の和集合をとる．
	@param &gene1 和する個体群１
	@param &gene2 和する個体群２
	@param &joinedGene 和集合した個体群
*/
void GaCommon::joinPopulation(
	const std::vector<std::vector<int> > &gene1,
	const std::vector<std::vector<int> > &gene2,
	std::vector<std::vector<int> > &joinedGene)
{
	// gene1の遺伝子をすべてコピー
	GaCommon::pushBackAll(joinedGene, gene1);

	for (auto match = gene2.begin(); match != gene2.end(); ++match)
	{
		auto obj	= std::find(joinedGene.begin(), joinedGene.end(), *match);
		if (obj == joinedGene.end())
			// gene2にのみ存在する遺伝子を追加
			joinedGene.push_back(*match);
	}
}

/*
	2次元配列の要素をもう一方の2次元配列にすべて格納する
	@param &targetPopulation 追加する対象の個体群
	@param &pushedPopulation 追加する要素を含む個体群
*/
void GaCommon::pushBackAll(
	std::vector<std::vector<int> > &targetPopulation,
	const std::vector<std::vector<int> > &pushedPopulation)
{
	targetPopulation.reserve(pushedPopulation.size());
	for (int tmp = 0; tmp < pushedPopulation.size(); ++tmp)
		targetPopulation.push_back(pushedPopulation[tmp]);
}