/***********
	GA.cpp

	規約（2016/05/09）
		- std::vector< std::vector<int> >の変数名は~Population
************/

#include "../include/Ga.h"
#include "../include/GaCommon.h"


GA::GA(int numVariable)
{
	std::cout << "Calls constructor." << std::endl;

	int geneLength	= 20;	// 遺伝子長
	int population		= 120;	// 個体数

	this->_population	= population;
	this->_geneLength	= geneLength;
	this->_numVariable	= numVariable;

	this->_searchPopulation	= std::vector< std::vector<int> >(population, std::vector<int>(geneLength));
}

GA::~GA()
{
}

/*
	初期化処理
*/
void GA::initGene()
{
	this->_initSearchPopulation();
	this->_archivePopulation	= std::vector< std::vector<int> >(this->_population, std::vector<int>(this->_geneLength*this->_numVariable));
}

/*
	探索母集団を初期化する．
*/
void GA::_initSearchPopulation()
{
	// カウント変数
	int tmp_column, tmp_row;

	int duplicateFlg;

	// [0.0, 1.0]のランダム値を作成
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<> randomValue(0.0, 1.0);

	do {
		duplicateFlg = 0;
		for (tmp_column = 0; tmp_column < this->_population; tmp_column++)
		{
			// 遺伝子の作成
			for (tmp_row = 0; tmp_row < this->_geneLength*this->_numVariable; tmp_row++)
			{
				this->_searchPopulation.at(tmp_column).at(tmp_row) = (randomValue(mt) > 0.5) ? 1 : 0;
			}

			if (this->_isDuplicatedGene(this->_searchPopulation, tmp_column))
				duplicateFlg = 1;
		}
	} while (duplicateFlg == 1);	
}

/*
	個体集団を表示する．
	テスト用
*/
void GA::_outputIndividuals(std::vector<std::vector<int>> &individuals)
{
	// カウント変数
	int tmp_column, tmp_row;

	for (tmp_column = 0; tmp_column < this->_population; tmp_column++)
	{
		for (tmp_row = 0; tmp_row < this->_geneLength*this->_numVariable; tmp_row++)
		{
			std::cout << individuals[tmp_column][tmp_row];
		}
		std::cout << std::endl;
	}
}

/*
	@param gene 個体集団の2次元配列
*/
bool GA::_isDuplicatedGene(std::vector<std::vector<int>> &gene, int column)
{
	// カウント変数
	int tmp_column, tmp_row;

	int duplicateFlg;
	bool result = false;

	for (tmp_column = 0; tmp_column < column; tmp_column++)
	{
		duplicateFlg = 1;
		for (tmp_row = 0; tmp_row < this->_geneLength; tmp_row++)
		{
			// 1つでも遺伝子が異なればフラグを回収
			if (gene[column][tmp_row] != gene[tmp_column][tmp_row])
				duplicateFlg = 0;
		}

		if (duplicateFlg == 1)
			result = true;
	}

	return result;
}

/*
	1個体の2進数データを表現型に変換
	@param &binary 1個体の遺伝子
	@param &phenotype 遺伝子型を表現型に変換したもの
*/
void GA::_convertPhenotype(const std::vector<int> &binary, std::vector<double> &phenotype)
{
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	for (int numVar = 0; numVar < this->_numVariable; ++numVar)
	{
		for (int numBinary = 0; numBinary < this->_geneLength; ++numBinary)
		{
			tmpGene.push_back(binary[numVar*this->_geneLength + numBinary]);
		}
		phenotype.push_back(this->_binary2Phenotype(tmpGene));
		tmpGene.clear();
	}
}

/*
	1遺伝子の2進数データを表現型に変換
	@param &binary 1遺伝子長の長さを持つ2進数データ
*/
double GA::_binary2Phenotype(const std::vector<int> &binary)
{
	// カウント変数
	int tmp;

	double decimal	= 0.;
	int place	= 0;

	for (tmp = this->_geneLength - 1; tmp >= 0 ; --tmp)
	{
		if (binary[tmp] == 1)
		{
			decimal	+= pow(2.0, (double)place);
			place	+= 1;
		}
	}

	return decimal / (pow(2.0, (double)this->_geneLength) - 1.0);
}

/*
	目的関数を計算
	@param &variable 変数
	@param &objectiveValue 目的関数の値（目的関数の数だけ領域を確保しておく）
*/
void GA::_getObjectiveFunc(const std::vector<double> &variable, std::vector<double> &objectiveValue)
{
	objectiveValue.at(0)	= this->_f1(variable);
	objectiveValue.at(1)	= this->_f2(variable);
}

double GA::_f1(const std::vector<double> &variable)
{
	return 1. - exp(-4.*variable[0])*pow(sin(6.*M_PI*variable[0]),6);
}

double GA::_f2(const std::vector<double> &var)
{
	double g, h;

	g = 1. + 9.*pow((var[1]+var[2]+var[3]+var[4]+var[5]+var[6]+var[7]+var[8]+var[9])/(9.),0.25);
	h = 1. - pow(this->_f1(var)/g,2);

	return g*h;
}

/*
	1個体から目的関数値を計算
	@param &binary 1個体の遺伝子
	@param &obj 目的関数の値（目的関数の数だけ領域を確保しておく）
*/
void GA::_binary2ObjectiveFunc(const std::vector<int> &binary, std::vector<double> &obj)
{
	// 表現型を一時的に保存
	std::vector<double> tmpPhenotype(this->_numVariable);

	this->_convertPhenotype(binary, tmpPhenotype);
	this->_getObjectiveFunc(tmpPhenotype, obj);
}

/***************
	未テスト関数
***************/

/*
	NSGA2実行用メソッド
*/
void GA::nsga2Run()
{
	int generation = 0;
	int maxGeneration	= 500;
	std::vector<std::vector<int> > margedPopulation, nextRankPopulation, crowdingSortedPopulation;
	std::vector<std::vector<std::vector<int> > > archivePopulation(maxGeneration), searchPopulation(maxGeneration), classifiedByRankGene;

	std::copy(this->_archivePopulation.begin(), this->_archivePopulation.end(), std::back_inserter(archivePopulation[0]));
	std::copy(this->_searchPopulation.begin(), this->_searchPopulation.end(), std::back_inserter(searchPopulation[0]));

	/*** Step1 ***/
	this->initGene();

	while(1)
	{
		/*** Step2 ***/
		// TODO:評価方法の確立
		this->_outputObjectiveValue(searchPopulation[generation], generation);

		/*** Step3 ***/
		GaCommon::joinPopulation(archivePopulation[generation], searchPopulation[generation], margedPopulation);
		this->_nonSuperioritySort(margedPopulation, classifiedByRankGene);
		margedPopulation.clear();

		/*** Step4 ***/
		this->_updateArchivePopulation(classifiedByRankGene, archivePopulation[generation+1], nextRankPopulation);

		/*** Step5 ***/
		this->_crowdingSort(nextRankPopulation, crowdingSortedPopulation);
		this->_insertIndividuals(archivePopulation[generation+1], nextRankPopulation);
		nextRankPopulation.clear();
		crowdingSortedPopulation.clear();

		/*** Step6 ***/
		if (generation < maxGeneration)
			break;

		/*** Step7 ***/
		// 選択と交叉を同時に行っている
		this->_crowdedTournamentSelection(archivePopulation[generation+1], searchPopulation[generation+1], classifiedByRankGene);
		classifiedByRankGene.clear();

		/*** Step8 ***/
		this->_mutationGene(searchPopulation[generation+1], 1/this->_geneLength*this->_numVariable);

		generation	+= 1;
	}
}

/*
	母集団に対して非優越ソートを行う
	@param gene ソートを行う母集団
	@param &classifiedByRankGene ランクごとにクラス分けした個体集団 
*/
void GA::_nonSuperioritySort(
	const std::vector <std::vector<int> > &targetPopulation,
	std::vector<std::vector<std::vector<int> > > &classifiedByRankGene)
{
	// カウント変数
	int tmp, tmp_column;

	int numGene, numObj, superiorityFlg;
	std::vector<double> targetObject(2), comparedObject(2);	// 目的関数の数
	std::vector<std::vector<int> > sortingPopulation, tmpRankedGene;

	// 個体にランクを付けていき，ランク付けされた個体は除く
	std::copy(targetPopulation.begin(), targetPopulation.end(), std::back_inserter(sortingPopulation));
	while(sortingPopulation.size() > 0)
	{
		numGene	= 0;
		while(numGene < sortingPopulation.size())
		{
			superiorityFlg	= 1;
			this->_binary2ObjectiveFunc(sortingPopulation[numGene], targetObject);

			for (tmp_column = 0; tmp_column < sortingPopulation.size(); ++tmp_column)
			{
				if (numGene == tmp_column)
					continue;
				this->_binary2ObjectiveFunc(sortingPopulation[tmp_column], comparedObject);

				for (numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
					if (targetObject[numObj] < comparedObject[numObj])
						superiorityFlg	= 0;
			}

			if (superiorityFlg == 1)
				tmpRankedGene.push_back(sortingPopulation[numGene]);
			targetObject.clear();
			numGene += 1;
		}

		// ランクごとに個体を保存し，個体群を更新
		classifiedByRankGene.push_back(tmpRankedGene);
		for (tmp = 0; tmp < tmpRankedGene.size(); ++tmp)
			GaCommonTemp<int>::removeElement(sortingPopulation, tmpRankedGene[tmp]);
		tmpRankedGene.clear();
	}
}

/*
	新たなアーカイブ集団を生成
	@param &classifiedByRankGene ランクごとの個体
	@param &newArchivePopulation ランク上位から取得した更新用アーカイブ母集団
	@param &nextRankPopulation アーカイブ母集団に入りきらなかった最高ランクの個体集団
*/
void GA::_updateArchivePopulation(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	std::vector<std::vector<int> > &newArchivePopulation,
	std::vector<std::vector<int> > &nextRankPopulation)
{
	int rank;

	for (rank = 0; rank < classifiedByRankGene.size(); ++rank)
	{
		if (this->_population - newArchivePopulation.size() > classifiedByRankGene[rank].size())
			GaCommon::pushBackAll(newArchivePopulation, classifiedByRankGene[rank]);
		else
			break;
	}

	// アーカイブ母集団に入りきらなかった最高ランクの個体集団を格納
	std::copy(classifiedByRankGene[rank].begin(), classifiedByRankGene[rank].end(), std::back_inserter(nextRankPopulation));
}

/*
	混雑度ソート
	@param &certainRankPopulation あるランクの個体群
	@param &crowdingSortedPopulation 混雑度ごとにソートされた個体群
*/
void GA::_crowdingSort(
	const std::vector<std::vector<int > > &certainRankPopulation,
	std::vector<std::vector<int> > &crowdingSortedPopulation)
{
	// カウント変数
	int tmp1, tmp2;

	int numObj;
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);
	std::vector<double> tmpObjectFunc(2);	// 目的関数の数
	std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// 目的関数の数	

	// 目的関数ごとに目的関数値が悪い順に並べる
	this->_putObjectiveSortedGeneEveryObjectiveFunc(certainRankPopulation, objectiveSortedGene);

	// 混雑度が大きい順に個体をソート
	std::copy(certainRankPopulation.begin(), certainRankPopulation.end(), std::back_inserter(crowdingSortedPopulation));
	double distance1 = 0., distance2 = 0.;
	for (tmp1 = 1; tmp1 < certainRankPopulation.size()-1; ++tmp1)
	{
		for (tmp2 = certainRankPopulation.size()-2;  tmp2 > tmp1; --tmp2)
		{
			distance1	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, crowdingSortedPopulation[tmp2]);
			distance2	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, crowdingSortedPopulation[tmp2-1]);
			if (distance1 > distance2)
			{
				tmpGene	= crowdingSortedPopulation[tmp2];
				crowdingSortedPopulation[tmp2]	= crowdingSortedPopulation[tmp2-1];
				crowdingSortedPopulation[tmp2-1]	= tmpGene;
			}
		}
	}
}

/*
	目的関数ごとに目的関数値が悪い順に並べる
	@param &targetPopulation ソートしたい個体群
	@param &objectiveSortedGene 目的関数ごとにソートされた個体（目的関数の数だけ領域を確保しておく）
*/
void GA::_putObjectiveSortedGeneEveryObjectiveFunc(
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<std::vector<std::vector<int> > > &objectiveSortedGene)
{
	int numObj;
	for (numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
	{
		this->_sortByObjectiveValue(targetPopulation, objectiveSortedGene[numObj], numObj);
		std::reverse(objectiveSortedGene[numObj].begin(), objectiveSortedGene[numObj].end());
	}
}

/*
	指定した目的関数値が小さい順に個体をバブルソートする
	@param &targetGene ソートしたい個体群
	@param &objectiveSortedGene ソート後の個体群
	@param num 対象とする目的関数の番号
*/
void GA::_sortByObjectiveValue(
	const std::vector<std::vector<int> > &targetGene,
	std::vector<std::vector<int> > &objectiveSortedPopulation,
	int num)
{
	// カウント変数
	int tmp1, tmp2;

	std::vector<double> tmpObject1(2);	// 目的関数の数
	std::vector<double> tmpObject2(2);	// 目的関数の数
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	std::copy(targetGene.begin(), targetGene.end(), std::back_inserter(objectiveSortedPopulation));

	for (tmp1 = 0; tmp1 < targetGene.size(); ++tmp1)
	{
		for (tmp2 = targetGene.size() - 1; tmp2 > tmp1; --tmp2)
		{
			this->_binary2ObjectiveFunc(targetGene[tmp2], tmpObject1);
			this->_binary2ObjectiveFunc(targetGene[tmp2-1], tmpObject2);
			if (tmpObject1[num] < tmpObject2[num])
			{
				tmpGene	= objectiveSortedPopulation[tmp2];
				objectiveSortedPopulation[tmp2]	= objectiveSortedPopulation[tmp2-1];
				objectiveSortedPopulation[tmp2-1]	= tmpGene;
			}
		}
	}
}

/*
	指定した目的関数に対して混雑度を計算する
	@param &objectiveSortedGene 目的関数ごとにソートされた個体群
	@param numGene 個体の番号（境界個体は選択できない）
	@param numObj 目的関数の番号
*/
double GA::_culcCrowdingDistanse(
	const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene,
	int numGene,
	int numObj)
{
	std::vector<double> tmpObjLeft(2);		// 目的関数の数
	std::vector<double> tmpObjRight(2);	// 目的関数の数
	std::vector<double> tmpObjMax(2);		// 目的関数の数
	std::vector<double> tmpObjMin(2);		// 目的関数の数

	// 境界個体を除いて個体を選択
	double distance	= 0.;
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][numGene-1], tmpObjLeft);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][numGene+1], tmpObjRight);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][0], tmpObjMax);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][objectiveSortedGene.front().size()-1], tmpObjMin);
	distance	= (tmpObjLeft[numObj] - tmpObjRight[numObj])/(tmpObjMax[numObj] - tmpObjMin[numObj]);

	return distance;
}

/*
	指定した個体の総混雑度を計算する
	@param &objectiveSortedGene 目的関数ごとに目的関数地が悪い順にソートされた個体群
	@param &individual 個体の遺伝子
*/
double GA::_culcCrowdingDistanseForIndividual(
	const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene,
	const std::vector<int> &individual)
{
	int numObj, numGene;
	std::vector<double> distance(2);	// 目的関数の数

	for (numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
	{
		for (numGene = 1; numGene < objectiveSortedGene.front().size()-1; ++numGene)
		{
			if (objectiveSortedGene[numObj][numGene] == individual)
			{
				distance[numObj]	= this->_culcCrowdingDistanse(objectiveSortedGene, numGene, numObj);
				break;
			}
		}
	}

	double resultDistanse	= 0.;
	for (numObj = 0; numObj < distance.size(); ++numObj)
	{
		resultDistanse	+= distance[numObj];
	}

	return  resultDistanse;
}

/*
	個体数がNになるまで個体を追加する
	@param &insertedPopulation 個体を追加される個体集団（個体数はN以下）
	@param &insertPopulation 追加する個体を含む個体集団
*/
void GA::_insertIndividuals(
	std::vector<std::vector<int> > &insertedPopulation,
	const std::vector<std::vector<int> > &insertPopulation)
{
	if (insertedPopulation.size() >= this->_population)
	{
		std::cout << "ERROR:Excess individual!!" << std::endl;
		return;
	}

	for (int tmp = 0; insertedPopulation.size() < this->_population; ++tmp)
		insertedPopulation.push_back(insertPopulation[tmp]);
}

/*
	混雑度トーナメント選択
	@param &selectedPopulation 選択される個体集団
	@param &newSearchPopulation 新たな探索母集団
	@param &classifiedByRankGene ランクごとにクラス分けされた個体
*/
void GA::_crowdedTournamentSelection(
	const std::vector<std::vector<int> > &selectedPopulation,
	std::vector<std::vector<int> > &newSearchPopulation,
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene)
{
	std::vector<int> parentGene1, parentGene2, childGene1, childGene2;
	std::vector<std::vector<int> > tmpSelectionPopulation, highRankPopulation;

	// 個体数がNになるまで選択を実行
	// 親個体をランダムに選択し，一様交叉を実行
	// 親×2と子×2のうちランクが上位の個体を1つ選択し，新たな探索母集団に追加する
	for (int numGene = 0; newSearchPopulation.size() < this->_population; ++numGene)
	{
		this->_select2GenesFromPopulation(selectedPopulation, parentGene1, parentGene2);
		this->_uniformCrossover(parentGene1, parentGene2, childGene1, childGene2);

		tmpSelectionPopulation.push_back(parentGene1);
		tmpSelectionPopulation.push_back(parentGene2);
		tmpSelectionPopulation.push_back(childGene1);
		tmpSelectionPopulation.push_back(childGene2);

		this->_highRankGeneSelection(classifiedByRankGene, tmpSelectionPopulation, highRankPopulation, 1);
		newSearchPopulation.push_back(highRankPopulation[1]);
	}
}

/*
	個体のランクを返す
	@param &classifiedByRankGene ランクごとにクラス分けされた個体
	@param &targetGene ランクを知りたい個体
*/
int GA::_returnGeneRank(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	const std::vector<int> &targetGene)
{
	int rank, numGene;

	for (rank = 0; rank < classifiedByRankGene.size(); ++rank)
	{
		for (numGene = 0; numGene < classifiedByRankGene[rank].size(); ++numGene)
		{
			if (classifiedByRankGene[rank][numGene] == targetGene)
				return rank;
		}
	}
}

/*
	個体群からランダムに2個体を選択する
	@param &targetPopulation 選択する個体群
	@param &gene1 選択された個体1
	@param &gene2 選択された個体2
*/
void GA::_select2GenesFromPopulation(
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<int> &gene1,
	std::vector<int> &gene2)
{
	int geneNum1, geneNum2;

	// 個体を選択するためのランダム値を生成
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, targetPopulation.size()-1);

	// ランダムに個体を選択
	do
	{
		geneNum1	= randomValue(mt);
		geneNum2	= randomValue(mt);
	} while (geneNum1 == geneNum2);

	gene1	= targetPopulation[geneNum1];
	gene2	= targetPopulation[geneNum2];
}

/*
	一様交叉をおこなう
	@param &parentGene1 親個体1
	@param &parentGene2 親個体2
	@param &childGene1 子個体1
	@param &childGene2 子個体2
*/
void GA::_uniformCrossover(
	const std::vector<int> &parentGene1,
	const std::vector<int> &parentGene2,
	std::vector<int> &childGene1,
	std::vector<int> &childGene2)
{
	// カウント変数
	int tmp;

	// マスクパターン用のランダム値生成器
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, 1);

	// マスクパターンを生成
	std::vector<int> maskPattern(this->_geneLength*this->_numVariable);
	for (tmp = 0; tmp < this->_geneLength; tmp++)
		maskPattern.push_back(randomValue(mt));

	// 交叉
	for (tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
	{
		if (maskPattern[tmp] == 0)
		{
			childGene1.push_back(parentGene1[tmp]);
			childGene2.push_back(parentGene2[tmp]);
		}
		else if (maskPattern[tmp] == 1)
		{
			childGene1.push_back(parentGene2[tmp]);
			childGene2.push_back(parentGene1[tmp]);
		}
	}
}

/*
	指定した個体から上位ランク個体を選択
	@param &classifiedByRankGene ランクごとにクラス分けされた個体
	@param &targetPopulation 選択する個体群
	@param &highRankPopulation 上位ランクの個体群
	@param num 選択する個体数
*/
void GA::_highRankGeneSelection(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<std::vector<int> > &highRankPopulation,
	int num)
{
	int geneRank;
	double geneDistance, longestDistance = 0.;
	std::vector<int> tmpHighRankGene(this->_geneLength*this->_numVariable);
	std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// 目的関数の数	

	// 対象の個体群をランクごとに分ける
	std::vector<std::vector<std::vector<int> > > tmpClassifiedByRankGene;
	this->_nonSuperioritySort(targetPopulation, tmpClassifiedByRankGene);

	// 上位ランクの個体をnum個選択する
	for (int rank = 0; highRankPopulation.size() < num; ++rank)
	{
		if (tmpClassifiedByRankGene[rank].size() == 1)
			highRankPopulation.push_back(tmpClassifiedByRankGene[rank][0]);
		else if (tmpClassifiedByRankGene.size() > 1)
		{
			// 同ランクに複数個体が存在した場合
			// 最も混雑距離が長い個体を選択して追加する
			for (int tmp = 0; tmp < tmpClassifiedByRankGene[rank].size(); ++tmp)
			{
				// 全個体中のランクから混雑度距離を求める
				geneRank	= this->_returnGeneRank(classifiedByRankGene, tmpClassifiedByRankGene[rank][tmp]);
				this->_putObjectiveSortedGeneEveryObjectiveFunc(classifiedByRankGene[geneRank], objectiveSortedGene);
				geneDistance	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, tmpClassifiedByRankGene[rank][tmp]);
				if (longestDistance < geneDistance)
				{
					longestDistance = geneDistance;
					tmpHighRankGene	= tmpClassifiedByRankGene[rank][tmp];
				}
			}

			highRankPopulation.push_back(tmpHighRankGene);
		}
	}
}

/*
	5番目の個体に突然変異を行う．
	適応度計算前の個体の遺伝子(0 or 1)をランダムに入れ替える
	@param &targetPopulation 対象の個体群
	@param mutationRate 突然変異率
*/
void GA::_mutationGene(
	std::vector<std::vector<int> > &targetPopulation,
	double mutationRate)
{
	int geneNum = 5;

	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);

	for (int tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
	{
		if (mutationRate > randomValue(mt))
		{
			switch (targetPopulation[geneNum][tmp])
			{
			case 0:
				targetPopulation[geneNum][tmp] = 1;
			case 1:
				targetPopulation[geneNum][tmp] = 0;
			default:
				break;
			}
		}
	}
}

/*
	指定した世代の個体集団のxと適応度を表示する
	@param &targetPopulation 対象の個体郡
	@param generation 現時点での世代
*/
void GA::_outputObjectiveValue(
	std::vector<std::vector<int> > targetPopulation,
	int generation)
{
	double x;
	std::vector<std::vector<double> > objectiveParameter(targetPopulation.size(), std::vector<double>(2));
	std::vector<double> objectiveValue(2);

	std::cout << generation << "-generation" << std::endl;
	for (int numGene = 0; numGene < targetPopulation.size(); ++numGene)
	{
		this->_binary2ObjectiveFunc(targetPopulation[numGene], objectiveParameter[numGene]);
		this->_getObjectiveFunc(objectiveParameter[numGene], objectiveValue);

		for (int numVar = 0; numVar < this->_numVariable; ++numVar)
		{
			std::cout << objectiveParameter[numGene][numVar] << ",";
		}
		std::cout << "\t";
		for (int numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
		{
			std::cout << objectiveValue[numObj] << ",";
		}
	}
	std::cout << std::endl;
}