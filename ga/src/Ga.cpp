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
	int population		= 50;	// 個体数

	this->_population	= population;
	this->_geneLength	= geneLength;
	this->_numVariable	= numVariable;
}

GA::~GA()
{
}

/*
	探索母集団を初期化する．
	@param &searchPopulation
*/
void GA::_initSearchPopulation(std::vector<std::vector<int> > &searchPopulation)
{
	while (searchPopulation.size() < this->_population)
	{
		std::vector<int> tmpGene(this->_geneLength*this->_numVariable);
		this->_createRandomlyIndividual(tmpGene);

		searchPopulation.reserve(this->_population);
		if (!this->_isDuplicatedGene(searchPopulation, tmpGene))
			searchPopulation.push_back(tmpGene);
	}
}

/*
	個体をランダムに1つ生成する
	@param &individual 生成する個体
*/
void GA::_createRandomlyIndividual(std::vector<int> &individual)
{
	// [0.0, 1.0]のランダム値を生成
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<> randomValue(0.0, 1.0);

	for (int numBinary = 0; numBinary < this->_geneLength*this->_numVariable; ++numBinary)
		individual[numBinary]	= (randomValue(mt) > 0.5) ? 1 : 0;
}

/*
	個体を表示する
	@param &individual 表示する個体
*/
void GA::_outputGene(const std::vector<int> &individual)
{
	int numBinary;

	for (numBinary = 0; numBinary < this->_geneLength*this->_numVariable; ++numBinary)
		std::cout << individual[numBinary];

	std::cout << std::endl;
}

/*
	個体集団を表示する
	@param &targetPopulation 表示する個体集団
*/
void GA::_outputPopulation(const std::vector<std::vector<int>> &targetPopulation)
{
	int numIndividual;

	for (numIndividual = 0; numIndividual < targetPopulation.size(); ++numIndividual)
		this->_outputGene(targetPopulation[numIndividual]);
}

/*
	個体が個体集団の中で重複していないか判定する
	@param &searchPopulation 個体集団
	@param &targetGene 個体
*/
bool GA::_isDuplicatedGene(std::vector<std::vector<int>> &searchPopulation, const std::vector<int> &targetGene)
{
	// 個体集団の要素が0だったら重複していない
	if (searchPopulation.size() == 0)
		return false;

	int numGene;
	for (numGene = 0; numGene < searchPopulation.size(); ++numGene)
		if (searchPopulation[numGene] == targetGene)
			return true;

	return false;
}

/*
	1個体2進数データを表現型に変換
	@param &binary 1個体の遺伝子
	@param &phenotype 遺伝子型を表現型に変換したもの（領域を確保しておく）
*/
void GA::_convertPhenotype(const std::vector<int> &binary, std::vector<double> &phenotype)
{

	for (int numVar = 0; numVar < this->_numVariable; ++numVar)
	{
		std::vector<int> tmpGene(this->_geneLength);

		for (int numBinary = 0; numBinary < this->_geneLength; ++numBinary)
			tmpGene[numBinary]	= binary[numVar*this->_geneLength + numBinary];

		phenotype[numVar]	= this->_binary2Phenotype(tmpGene);
	}
}

/*
	1遺伝子の2進数データを表現型に変換
	@param &binary 1遺伝子長の長さを持つ2進数データ
*/
double GA::_binary2Phenotype(const std::vector<int> &binary)
{
	int place		= 0;
	double decimal	= 0.;

	for (int numBinary = this->_geneLength - 1; numBinary >= 0 ; --numBinary)
	{
		if (binary[numBinary] == 1)
			decimal	+= pow(2.0, place);
		place	+= 1;
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
	objectiveValue[0]	= this->_f1(variable);
	objectiveValue[1]	= this->_f2(variable);
}

double GA::_f1(const std::vector<double> &variable)
{
	return 1. - exp(-4.*variable[0])*pow(sin(6.*M_PI*variable[0]),6);
	// return -(sin(3.*variable[0])+0.5*sin(9.*variable[0])+sin(15.*variable[0]+50.));
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
//	GaCommonTemp<double>::outputAllElement(obj);
}

/*
	NSGA2実行用メソッド
*/
void GA::nsga2Run()
{
	int generation = 0;
	int maxGeneration	= 500;
	std::vector<std::vector<std::vector<int> > > archivePopulation(maxGeneration), searchPopulation(maxGeneration);

	/*** Step1 ***/
	std::cout << "Call Step1!" << std::endl;
	this->_initSearchPopulation(searchPopulation[generation]);

	while(1)
	{
		/*** Step2 ***/
		std::cout << "Call Step2!" << std::endl;
		// TODO:評価方法の確立
//		this->_outputObjectiveValue(searchPopulation[generation], generation);

		/*** Step3 ***/
		std::cout << "Call Step3!" << std::endl;
		std::vector<std::vector<std::vector<int> > > classifiedByRankGene;
		std::vector<std::vector<int> > margedPopulation;
		GaCommon::joinPopulation(archivePopulation[generation], searchPopulation[generation], margedPopulation);
//		std::cout << margedPopulation.size() << std::endl;
		this->_nonSuperioritySort(margedPopulation, classifiedByRankGene);

//		std::cout << classifiedByRankGene[0].size() << std::endl;

		/*** Step4 ***/
		std::cout << "Call Step4!" << std::endl;
		std::vector<std::vector<int> > nextRankPopulation;
		this->_updateArchivePopulation(classifiedByRankGene, archivePopulation[generation+1], nextRankPopulation);
//		std::cout << archivePopulation[generation+1].size() << std::endl;

		/*** Step5 ***/
		std::cout << "Call Step5!" << std::endl;
		std::vector<std::vector<int> > crowdingSortedPopulation;
		this->_crowdingSort(nextRankPopulation, crowdingSortedPopulation);
		this->_insertIndividuals(archivePopulation[generation+1], nextRankPopulation);
//		std::cout << archivePopulation[generation+1].size() << std::endl;

		/*** Step6 ***/
		std::cout << "Call Step6!" << std::endl;
		this->_outputObjectiveValue(searchPopulation[generation], generation);
		if  (generation == maxGeneration)
			break;

		/*** Step7 ***/
		std::cout << "Call Step7!" << std::endl;
		// 選択と交叉を同時に行っている
		this->_crowdedTournamentSelection(archivePopulation[generation+1], searchPopulation[generation+1], classifiedByRankGene);
//		std::cout << searchPopulation[generation+1].size() << std::endl;

		/*** Step8 ***/
		std::cout << "Call Step8!" << std::endl << std::endl;
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
	int tmp, numGene;
	std::vector<std::vector<int> > sortingPopulation;

	// 個体にランクを付けていき，ランク付けされた個体は除く
	std::copy(targetPopulation.begin(), targetPopulation.end(), std::back_inserter(sortingPopulation));
	while (sortingPopulation.size() > 0)
	{
		std::vector<std::vector<int> > tmpRankedGene;

		for (numGene = 0; numGene < sortingPopulation.size(); ++numGene)
		{
//			std::cout << this->_isSuperior(sortingPopulation[numGene], sortingPopulation) << std::endl;
			if (this->_isSuperior(sortingPopulation[numGene], sortingPopulation))
				tmpRankedGene.push_back(sortingPopulation[numGene]);			
		}

		// ランクごとに個体を保存し，個体群を更新
		classifiedByRankGene.push_back(tmpRankedGene);
		for (tmp = 0; tmp < tmpRankedGene.size(); ++tmp)
			GaCommonTemp<int>::removeElement(sortingPopulation, tmpRankedGene[tmp]);
	}
}

/*
	個体が個体集団の中で優越しているかを判定する
	@param &targetGene 判定する個体
	@param &comparedPopulation targetGeneが属する個体集団
*/
bool GA::_isSuperior(
	const std::vector<int> &targetGene,
	const std::vector<std::vector<int> > &comparedPopulation)
{
	int numGene, numSuperior;
	bool isSuperior;
	std::vector<double> targetObjeciveValue(2), comparedObjectiveValue(2);	// 目的関数の数

	isSuperior	= true;
	this->_binary2ObjectiveFunc(targetGene, targetObjeciveValue);
//	GaCommonTemp<double>::outputAllElement(targetObjeciveValue);

	for (numGene = 0; numGene < comparedPopulation.size(); ++numGene)
	{
		if (comparedPopulation[numGene] == targetGene)
			continue;

		this->_binary2ObjectiveFunc(comparedPopulation[numGene], comparedObjectiveValue);

		// 目的関数値を比較して1つも小さくなかったら非優越固体でない
		numSuperior	= this->_numOfSuperior(targetObjeciveValue, comparedObjectiveValue);
		if (numSuperior == 0)
			isSuperior	= false;
	}

	return isSuperior;
}

/*
	目的関数のうち，優越している個数を返す
	@param &targetObjectiveValue 比べる固体の目的関数値
	@param &comparedObjectiveValue 比べられる固体の目的関数値
*/
int GA::_numOfSuperior(
	const std::vector<double> &targetObjeciveValue,
	const std::vector<double> &comparedObjectiveValue)
{
	int num = 0;

	for (int numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
		if (targetObjeciveValue[numObj] < comparedObjectiveValue[numObj])
			num	+= 1;

	return num;
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
//		std::cout << classifiedByRankGene[rank].size() << std::endl;
		if (this->_population - newArchivePopulation.size() > classifiedByRankGene[rank].size())
			GaCommon::pushBackAll(newArchivePopulation, classifiedByRankGene[rank]);
		else
			break;
	}

	// アーカイブ母集団に入りきらなかった最高ランクの個体集団を格納
	nextRankPopulation.reserve(classifiedByRankGene[rank].size());
	nextRankPopulation	= classifiedByRankGene[rank];
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
	std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// 目的関数の数

	// 目的関数ごとに目的関数値が悪い順に並べる
	this->_putObjectiveSortedGeneEveryObjectiveFunc(certainRankPopulation, objectiveSortedGene);

	// 混雑度が大きい順に個体をソート
	crowdingSortedPopulation.reserve(certainRankPopulation.size());
	crowdingSortedPopulation	= certainRankPopulation;
	this->_orderByCrowdingDistanceUsingBubbleSort(crowdingSortedPopulation, objectiveSortedGene);
}

/*
	混雑度が大きい順に個体群をバブルソート
	@param &sortingPopulation ソート対象の個体群
	@param &objectiveSortedGene 目的関数ごとに目的関数値が悪い順に並んだ個体群
*/
void GA::_orderByCrowdingDistanceUsingBubbleSort(
	std::vector<std::vector<int> > &sortingPopulation,
	const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene)
{
	double distance1, distance2;
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	for (int tmp1 = 1; tmp1 < sortingPopulation.size(); ++tmp1)
	{
		for (int tmp2 = sortingPopulation.size()-1; tmp2 > tmp1 ; --tmp2)
		{
			distance1	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, sortingPopulation[tmp2]);
			distance2	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, sortingPopulation[tmp2-1]);
			if (distance1 > distance2)
			{
				tmpGene	= sortingPopulation[tmp2];
				sortingPopulation[tmp2]	= sortingPopulation[tmp2-1];
				sortingPopulation[tmp2-1]	= tmpGene;
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
	for (int numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
	{
		this->_orderBySmallObjectiveValueUsingBubbleSort(targetPopulation, objectiveSortedGene[numObj], numObj);
		std::reverse(objectiveSortedGene[numObj].begin(), objectiveSortedGene[numObj].end());
	}
}

/*
	指定した目的関数値が小さい順に個体をバブルソートする
	@param &targetGene ソートしたい個体群
	@param &objectiveSortedGene ソート後の個体群
	@param num 対象とする目的関数の番号
*/
void GA::_orderBySmallObjectiveValueUsingBubbleSort(
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<std::vector<int> > &objectiveSortedPopulation,
	int num)
{
	std::vector<double> tmpObject1(2), tmpObject2(2);	// 目的関数の数
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	// ソートする個体群をコピー
	objectiveSortedPopulation.reserve(targetPopulation.size());
	objectiveSortedPopulation	= targetPopulation;

	for (int tmp1 = 0; tmp1 < targetPopulation.size(); ++tmp1)
	{
		for (int tmp2 = targetPopulation.size() - 1; tmp2 > tmp1; --tmp2)
		{
			this->_binary2ObjectiveFunc(targetPopulation[tmp2], tmpObject1);
			this->_binary2ObjectiveFunc(targetPopulation[tmp2-1], tmpObject2);
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
	std::vector<double> tmpObjRight(2);		// 目的関数の数
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
		resultDistanse	+= distance[numObj];

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

	insertedPopulation.reserve(this->_population);
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
	std::vector<int> parentGene1(this->_geneLength*this->_numVariable), parentGene2(this->_geneLength*this->_numVariable);
	std::vector<int> childGene1(this->_geneLength*this->_numVariable), childGene2(this->_geneLength*this->_numVariable);

	// 個体数がNになるまで選択を実行
	// 親個体をランダムに選択し，一様交叉を実行
	// 親×2と子×2のうちランクが上位の個体を1つ選択し，新たな探索母集団に追加する
	newSearchPopulation.reserve(this->_population);
	for (int numGene = 0; newSearchPopulation.size() < this->_population; ++numGene)
	{
		std::vector<std::vector<int> > tmpSelectionPopulation(4), highRankPopulation;

		this->_select2GenesFromPopulation(selectedPopulation, parentGene1, parentGene2);
		this->_uniformCrossover(parentGene1, parentGene2, childGene1, childGene2);

		tmpSelectionPopulation[0]	= parentGene1;
		tmpSelectionPopulation[1]	= parentGene2;
		tmpSelectionPopulation[2]	= childGene1;
		tmpSelectionPopulation[3]	= childGene2;

		this->_highRankGeneSelection(classifiedByRankGene, tmpSelectionPopulation, highRankPopulation, 1);
		newSearchPopulation.push_back(highRankPopulation[0]);
//		std::cout << newSearchPopulation.size() << std::endl;
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
	@param &childGene1 子個体1（領域確保済み）
	@param &childGene2 子個体2（領域確保済み）
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
	for (tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
		maskPattern[tmp]	= randomValue(mt);

	// 交叉
	for (tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
	{
		if (maskPattern[tmp] == 0)
		{
			childGene1[tmp]	= parentGene1[tmp];
			childGene2[tmp]	= parentGene2[tmp];
		}
		else if (maskPattern[tmp] == 1)
		{
			childGene1[tmp]	= parentGene2[tmp];
			childGene2[tmp]	= parentGene1[tmp];
		}
	}
}

/*
	2点交叉を行う
	@param &parentGene1 親個体1
	@param &parentGene2 親個体2
	@param &childGene1 子個体1（領域確保済み）
	@param &childGene2 子個体2（領域確保済み）
*/
void GA::_2pointCrossover(
	const std::vector<int> &parentGene1,
	const std::vector<int> &parentGene2,
	std::vector<int> &childGene1,
	std::vector<int> &childGene2)
{
	// ランダム値生成器
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, this->_geneLength-1);

	int point1, point2;
	do
	{
		point1	= randomValue(mt);
		point2	= randomValue(mt);
	} while (point1 == point2);

	// 交叉点を決める
	int start, end;
	if (point1 < point2)
	{
		start	= point1;
		end	= point2;
	} else if (point1 > point2)
	{
		start	= point2;
		end	= point1;
	} else
	{
		std::cout << "ERROR:No crossover point!";
		return;
	}

	childGene1	= parentGene1;
	childGene2	= parentGene2;
	// this->_copyGene(parentGene1, childGene1);
	// this->_copyGene(parentGene2, childGene2);
	for (int tmp = start; tmp <= end; ++tmp)
	{
		childGene1[tmp]	= parentGene2[tmp];
		childGene2[tmp]	= parentGene1[tmp];
	}
}

/*
	遺伝子をコピーする
	@param &targetGene コピーする遺伝子
	@param &copiedGene コピーされる遺伝子（領域確保済み）
*/
void GA::_copyGene(
	const std::vector<int> &targetGene,
	std::vector<int> &copiedGene)
{
	for (int tmp = 0; tmp < targetGene.size(); ++tmp)
		copiedGene[tmp]	= targetGene[tmp];
}

/*
	指定した個体群から上位ランク個体を選択
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

	// 対象の個体群をランクごとに分ける
	std::vector<std::vector<std::vector<int> > > tmpClassifiedByRankGene;
	this->_nonSuperioritySort(targetPopulation, tmpClassifiedByRankGene);

	// 上位ランクの個体をnum個選択する
	highRankPopulation.reserve(num);
	for (int rank = 0; highRankPopulation.size() < num; ++rank)
	{
		if (tmpClassifiedByRankGene[rank].size() == 1)
			highRankPopulation.push_back(tmpClassifiedByRankGene[rank].front());
		else if (tmpClassifiedByRankGene[rank].size() > 1)
		{
			// 同ランクに複数個体が存在した場合
			// 最も混雑距離が長い個体を選択して追加する
			// TODO:ランクが高いほうが優先されるから複数あったら混雑距離の順にすべての個体を格納
			for (int tmp = 0; tmp < tmpClassifiedByRankGene[rank].size(); ++tmp)
			{
				// 全個体中のランクから混雑度距離を求める
				std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// 目的関数の数
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
	指定した世代の個体集団の目的関数値を表示する
	@param &targetPopulation 対象の個体郡
	@param generation 現時点での世代
*/
void GA::_outputObjectiveValue(
	std::vector<std::vector<int> > targetPopulation,
	int generation)
{
	double x;
	std::vector<double> objectiveValue(2);

	double tmp;
	std::cout << generation << "-generation" << std::endl;
	for (int numGene = 0; numGene < targetPopulation.size(); ++numGene)
	{
		this->_binary2ObjectiveFunc(targetPopulation[numGene], objectiveValue);

		// tmp	= 1. - pow(objectiveValue[0],2);
		// std::cout << objectiveValue[0] << " " << tmp << std::endl;

		for (int numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
		{
			std::cout << objectiveValue[numObj];
			if (numObj != 1)	// 目的関数の数 - 1
				std::cout << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/*
	指定した世代の個体集団のxを表示する
	@param &targetPopulation 対象の個体群
	@param generation 現時点での世代
*/
// void GA::_outputParameterValue(
// 	const std::vector<std::vector<int> > &targetPopulation,
// 	int generation)
// {
// 	for (int numGene = 0; numGene < targetPopulation.size(); ++numGene)
// 	{
// 		this->_binary2ObjectiveFunc
// 	}
// }