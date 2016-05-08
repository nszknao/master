#include "GaCommon.h"

/*
	２つの二次元配列の和集合をとる．
	３つ目の引数に結果を格納．
*/
void GaCommon::joinGene(const std::vector<std::vector<int> > &gene1, const std::vector<std::vector<int> > &gene2, std::vector<std::vector<int> > &joinedGene)
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

