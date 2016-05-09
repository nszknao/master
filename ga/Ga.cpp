/***********
	GA.cpp

	規約（2016/05/09）
		- std::vector< std::vector<int> >の変数名は~Population
************/

#include "Ga.h"
#include "GaCommon.h"


GA::GA(int numVariable)
{
	std::cout << "Calls constructor." << std::endl;

	int geneLength	= 20;	// 遺伝子長
	int population	= 120;	// 個体数

	this->_population	= population;
	this->_geneLength	= geneLength;
	this->_numVariable	= numVariable;

	this->_searchPopulation	= std::vector< std::vector<int> >(population, std::vector<int>(geneLength));
	this->fitness			= std::vector<int>(population);
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
	2つ目の引数に結果を格納
	@param &binary 1個体の遺伝子
	@param &phenotype 遺伝子型を表現型に変換したもの
*/
void GA::_convertPhenotype(const std::vector<int> &binary, std::vector<double> &phenotype)
{
	// カウント変数
	int tmp_var, tmp_cnt;

	std::vector<int> tmpGene(this->geneLength);

	for (tmp_cnt = 0; tmp_cnt < this->_numVariable; ++tmp_cnt)
	{
		for (tmp_var = 0; tmp_var < this->_geneLength; ++tmp_var)
		{
			tmpGene.push_back(binary[tmp_column][tmp_cnt*this->_geneLength + tmp_var]);
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
	@param &var 変数
	@param &obj 目的関数の値（目的関数の数だけ領域を確保しておく）
*/
void GA::_getObjectiveFunc(const std::vector<double> &var, std::vector<double> &obj)
{
	obj.at(0)	= this->_f1(var);
	obj.at(1)	= this->_f2(var);
}

double GA::_f1(const std::vector<double> &var)
{
	return 1. - exp(-4.*var[0])*pow(sin(6.*M_PI*var[0]),6);
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
	/*** Step1 ***/
	this->_initGene();

	/*** Step2 ***/
	// TODO:評価方法の確立

	/*** Step3 ***/
	std::vector<std::vector<int> > margedPopulation;
	GaCommon::joinPopulation(this->_archivePopulation, this->_searchPopulation, margedPopulation);
	std::vector<std::vector<std::vector<int> > > classifiedByRankGene;
	this->_nonSuperioritySort(margedPopulation, classifiedByRankGene);

	/*** Step4 ***/
	std::vector<std::vector<int> > newArchivePopulation, nextRankPopulation;
	this->_updateArchivePopulation(classifiedByRankGene, newArchivePopulation, nextRankPopulation);

	/*** Step5 ***/
	std::vector<std::vector<int> > crowdingSortedPopulation;
	this->_crowdingSort(nextRankPopulation, crowdingSortedPopulation);
	this->_insertIndividuals(newArchivePopulation, nextRankPopulation);

	/*** Step6 ***/
	// TODO:終了条件を設定

	/*** Step7 ***/
	
	/*** Step8 ***/
}

/*
	母集団に対して非優越ソートを行う
	@param gene ソートを行う母集団
	@param &classifiedByRankGene ランクごとにクラス分けした個体集団 
*/
void GA::_nonSuperioritySort(const std::vector <std::vector<int> > &targetPopulation, std::vector<std::vector<std::vector<int> > > &classifiedByRankGene)
{
	// カウント変数
	int tmp, tmp_column;

	int num, numObj, superiorityFlg;
	std::vector<double> targetObject(2), comparedObject(2);	// 目的関数の数
	std::vector<std::vector<int> > sortingPopulation, tmpRankedGene;

	// 個体にランクを付けていき，ランク付けされた個体は除く
	std::copy(targetPopulation.begin(), targetPopulation.end(), std::back_inserter(sortingPopulation));

	while(sortingPopulation.size() > 0)
	{
		num	= 0;
		while(num < sortingPopulation.size())
		{
			superiorityFlg	= 1;
			this->_binary2ObjectiveFunc(sortingPopulation[num], targetObject);

			for (tmp_column = 0; tmp_column < objValue.size(); ++tmp_column)
			{
				if (num == tmp_column)
					continue;
				this->_binary2ObjectiveFunc(sortingPopulation[tmp_column], comparedObject);

				for (numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
					if (targetObject[numObj] < comparedObject[numObj])
						superiorityFlg	= 0;
			}

			if (superiorityFlg == 1)
				tmpRankedGene.push_back(sortingPopulation[num]);
			targetObject.clear();
			num += 1;
		}

		// ランクごとに個体を保存し，個体群を更新
		classifiedByRankGene.push_back(tmpRankedGene);
		for (tmp = 0; tmp < tmpRankedGene.size(); ++tmp)
			GaCommon::removeElement(sortingPopulation, tmpRankedGene[tmp]);
		tmpRankedGene.clear();
	}
}

/*
	新たなアーカイブ集団を生成
	@param &classifiedByRankGene ランクごとの個体
	@param &newArchivePopulation ランク上位から取得した更新用アーカイブ母集団
	@param &nextRankPopulation アーカイブ母集団に入りきらなかった最高ランクの個体集団
*/
void GA::_updateArchivePopulation(const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene, std::vector<std::vector<int> > &newArchivePopulation, std::vector<std::vector<int> > &nextRankPopulation)
{
	for (int rank = 0; rank < classifiedByRankGene.size(); ++rank)
	{
		if (this->_population - newArchivePopulation.size() > classifiedByRankGene[rank].size())
			newArchivePopulation.push_back(classifiedByRankGene[rank]);
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
void GA::_crowdingSort(const std::vector<std::vector<int > > &certainRankPopulation, std::vector<std::vector<int> > &crowdingSortedPopulation)
{
	// カウント変数
	int tmp1, tmp2;

	int numObj;
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);
	std::vector<double> tmpObjectFunc(2);	// 目的関数の数
	std::vector<std::vector<std::vector<int> > > sortedObjectiveFunc(2);	// 目的関数の数	

	// 目的関数ごとに目的関数値が悪い順に並べる
	for (numObj = 0; numObj < 2; ++numObj)	// 目的関数の数
	{
		this->_sortRegardingObjective(certainRankPopulation, sortedObjectiveFunc[numObj], numObj);
		std::reverse(sortedObjectiveFunc[numObj].begin(), sortedObjectiveFunc[numObj].end());
	}

	// 混雑度が大きい順に固体をソート
	std::copy(certainRankPopulation.begin(), certainRankPopulation.end(), std::back_inserter(crowdingSortedPopulation));
	double distance1 = 0., distance2 = 0.;
	for (tmp1 = 1; tmp1 < certainRankPopulation.size()-1; ++tmp1)
	{
		for (tmp2 = certainRankPopulation.size()-2;  tmp2 > tmp1; --tmp2)
		{
			distance1	= this->_culcCrowdingDistanseForIndividual(sortedObjectiveFunc, crowdingSortedPopulation[tmp2]);
			distance2	= this->_culcCrowdingDistanseForIndividual(sortedObjectiveFunc, crowdingSortedPopulation[tmp2-1]);
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
	指定した目的関数値が小さい順に個体をバブルソートする
	@param &targetGene ソートしたい個体群
	@param &sortedGene ソート後の個体群
	@param num 対象とする目的関数の番号
*/
void GA::_sortRegardingObjective(const std::vector<std::vector<int> > &targetGene, std::vector<std::vector<int> > &sortedGene, int num)
{
	// カウント変数
	int tmp1, tmp2;

	std::vector<double> tmpObject1(2);	// 目的関数の数
	std::vector<double> tmpObject2(2);	// 目的関数の数
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	std::copy(targetGene.begin(), targetGene.end(), std::back_inserter(sortedGene));

	for (tmp1 = 0; tmp1 < targetGene.size(); ++tmp1)
	{
		for (tmp2 = targetGene.size() - 1; tmp2 > tmp1; --tmp2)
		{
			this->_binary2ObjectiveFunc(targetGene[tmp2], tmpObject1);
			this->_binary2ObjectiveFunc(targetGene[tmp2-1], tmpObject2);
			if (tmpObject1[num] < tmpObject2[num])
			{
				tmpGene	= sortedGene[tmp2];
				sortedGene[tmp2]	= sortedGene[tmp2-1];
				sortedGene[tmp2-1]	= tmpGene;
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
double GA::_culcCrowdingDistanse(const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene, int numGene, int numObj)
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
	@param &objectiveSortedGene 目的関数ごとにソートされた個体群
	@param &individual 個体の遺伝子
*/
double GA::_culcCrowdingDistanseForIndividual(const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene, const std::vector<int> &individual)
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
void GA::_insertIndividuals(std::vector<std::vector<int> > &insertedPopulation, const std::vector<std::vector<int> > &insertPopulation)
{
	if (insertedPopulation.size() >= this->_population)
	{
		std::cout << "ERROR:Excess individual!!" << std::endl;
		return;
	}

	for (int tmp = 0; insertedPopulation.size() < this->_population; ++tmp)
	{
		insertedPopulation.push_back(insertPopulation[tmp]);
	}
}




/***************
	古い関数
***************/

/*
	指定した世代の個体集団のxと適応度を表示する。
*/
void GA::outputGeneration(int generation)
{
	size_t tmp_column;

	double x;

	std::cout << generation << "-generation" << std::endl;
	std::cout << "number\tx\tfitness" << std::endl;
	for (tmp_column = 0; tmp_column < this->_geneLength*this->_numVariable; tmp_column++)
	{
		x = this->_binary2Phenotype(allIndividual[tmp_column]);
		std::cout << tmp_column << "\t" << x << "\t" << this->_getObjectiveFunc(x) << std::endl;
	}
}

void GA::uniformCrossover()
{
	// カウント変数
	int tmp;

	int parent1, parent2;
	std::vector<int> maskPattern, parent1Individual, parent2Individual, child1Individual, child2Individual;

	// 親を選択
	do
	{
		parent1 = this->selectIndividual();
		parent2 = this->selectIndividual();
		if (parent1 == -1 || parent2 == -1)
		{
			std::cout << "Couldn't select individual." << std::endl;
			return;
		}
	} while (parent1 == parent2);

	return;
	
	for (tmp = 0; tmp < this->_geneLength; tmp++)
	{
		parent1Individual.push_back(this->allIndividual[parent1][tmp]);
		parent2Individual.push_back(this->allIndividual[parent2][tmp]);
	}

	// マスクパターンを生成
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, 1);
	for (tmp = 0; tmp < this->_geneLength; tmp++)
	{
		maskPattern.push_back(randomValue(mt));
	}

	// 交叉
	for (tmp = 0; tmp < this->_geneLength; tmp++)
	{
		if (maskPattern[tmp] == 0)
		{
			child1Individual.push_back(parent1Individual[tmp]);
			child2Individual.push_back(parent2Individual[tmp]);
		}
		else if (maskPattern[tmp] == 1)
		{
			child1Individual.push_back(parent2Individual[tmp]);
			child2Individual.push_back(parent1Individual[tmp]);
		}
	}

	// 世代交代
	this->allIndividual[parent1].swap(child1Individual);
	this->allIndividual[parent2].swap(child2Individual);
}

/*
	適応度比例戦略から固体を選択し，その番号を返す．
	失敗した場合は-1を返す．
*/
int GA::selectIndividual()
{
	// カウント変数
	int tmp_column;

	std::vector<double> shuffledFitness, fitnessRatio, selectedIndividual;

	// 固体集団をシャッフルして順番を変えずに適応度を求める
	double sumFitness = 0.0, x;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(this->allIndividual.begin(), this->allIndividual.end(), std::mt19937(seed));
	for (tmp_column = 0; tmp_column < this->_population; tmp_column++)
	{
		x = this->_binary2Phenotype(this->allIndividual[tmp_column]);
		shuffledFitness.push_back(this->_getObjectiveFunc(x));
		sumFitness += this->_getObjectiveFunc(x);
	}

	// 適応度の比
	for (tmp_column = 0; tmp_column < this->_population; tmp_column++)
	{
		fitnessRatio.push_back(shuffledFitness[tmp_column] / sumFitness);
	}

	// 適応度から確率的に固体を選択
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);
	int selectedFlg = 0, individualNum = 0;
	do
	{
		for (individualNum = 0; individualNum < fitnessRatio.size(); individualNum++)
		{
			if (fitnessRatio[individualNum] > randomValue(mt))	// 毎回ランダム値が生成される	
			{
				return individualNum;
			}
		}
	} while (selectedFlg == 0);

	return -1;	// 固体が存在しなかった場合
}

/*
	ランキング法をもとに選択淘汰を行う．
*/
void GA::selectRanking()
{
	// カウント変数
	int tmp, tmp1, tmp2;

	// ソート順を記憶するためのキー
	std::vector<int> key;
	key.resize(this->_population);
	for (tmp = 0; tmp < this->_population; tmp++)
		key[tmp] = tmp;

	// 適応度を降順にソート．その順番をkeyに保存．
	double tmpVar;
	int tmpKey;
	for (tmp1 = 0; tmp1 < this->fitness.size(); ++tmp1)
	{
		for (tmp2 = tmp1+1; tmp2 < this->fitness.size(); ++tmp2)
		{
			if (this->fitness[tmp1] < this->fitness[tmp2])
			{
				tmpVar = this->fitness[tmp1];
				this->fitness[tmp1] = this->fitness[tmp2];
				this->fitness[tmp2] = tmpVar;

				tmpKey		= key[tmp1];
				key[tmp1]	= key[tmp2];
				key[tmp2]	= tmpKey;
			}
		}
	}

	// 個体集団を適応度の降順に並べる
	for (tmp = 0; tmp < this->_geneLength; tmp++)
	{
		this->allIndividual[tmp] = this->allIndividual[key[tmp]];
	}

	// 適応度が上位2つの個体をコピーし，下位2つを淘汰する．
	this->allIndividual[this->_geneLength - 1]	= this->allIndividual[0];
	this->allIndividual[this->_geneLength - 2]	= this->allIndividual[0];
}

/*
	5番目の個体に突然変異を行う．
	適応度計算前の個体の遺伝子(0 or 1)をランダムに入れ替える
	@param mutationRate 突然変異率
*/
void GA::mutation(double mutationRate)
{
	// カウント変数
	size_t tmp_row;

	int column = 5;

	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);

	for (tmp_row = 0; tmp_row < this->_geneLength; tmp_row++)
	{
		if (mutationRate > randomValue(mt))
		{
			switch (this->allIndividual[column][tmp_row])
			{
			case 0:
				this->allIndividual[column][tmp_row] = 1;
			case 1:
				this->allIndividual[column][tmp_row] = 0;
			default:
				break;
			}
		}
	}
}
