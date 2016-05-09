#include "GaCommon.h"

/*
	２つの二次元配列の和集合をとる．
	@param &gene1 和する個体群１
	@param &gene2 和する個体群２
	@param &joinedGene 和集合した個体群
*/
void GaCommon::joinPopulation(const std::vector<std::vector<int> > &gene1, const std::vector<std::vector<int> > &gene2, std::vector<std::vector<int> > &joinedGene)
{
	// gene1の遺伝子をすべてコピー
	std::copy(gene1.begin(), gene1.end(), std::back_inserter(joinedGene));

	for (auto match = gene2.begin(); match != gene2.end(); ++match)
	{
		auto obj	= std::find(joinedGene.begin(), joinedGene.end(), *match);
		if (obj == joinedGene.end())
		{
			// gene2にのみ存在する遺伝子を追加
			joinedGene.push_back(*match);
		}
	}
}

