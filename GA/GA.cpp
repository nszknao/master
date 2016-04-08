#include "GA.h"


GA::GA()
{
}

GA::~GA()
{
}

void GA::initGene()
{
	// カウント変数
	size_t tmp, tmp_column, tmp_row;

	int duplicateFlg;

	// [0.0, 1.0]のランダム値を作成
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);

	do {
		duplicateFlg = 0;
		for (tmp_column = 0; tmp_column < this->population; tmp_column++)
		{
			// 遺伝子の作成
			for (tmp_row = 0; tmp_row < this->geneLength; tmp_row++)
			{
				allIndividual[tmp_column][tmp_row] = (randomValue(mt) > 0.5) ? 1 : 0;
			}

			if (this->_isDuplicatedGene(allIndividual, tmp_column))
				duplicateFlg = 1;
		}
	} while (duplicateFlg = 1)
}

/*
	@param gene 個体集団の2次元配列
*/
bool GA::_isDuplicatedGene(std::vector<std::vector<int>> gene, int column)
{
	// カウント変数
	size_t tmp_column, tmp_row;

	int duplicateFlg;
	bool result = false;

	for (tmp_column = 0; tmp_column < column; tmp_column++)
	{
		duplicateFlg = 1;
		for (tmp_row = 0; tmp_row < this->geneLength; tmp_row++)
		{
			// 1つでも遺伝子が異なればフラグを回収
			if (gene[column][tmp_row] != gene[tmp_column][tmp_row])
				duplicateFlg = 0;
		}

		if (duplicateFlg = 1)
			result = true;
	}

	return result;
}

void GA::culcFitness()
{
	// カウント変数
	int tmp_column, tmp_row;

	double x;
	for (tmp_column = 0; tmp_column < this->population; tmp_column++)
	{
		x = this->_binary2Phenotype(allIndividual[tmp_column]);
		this->fitness[tmp_column] = this->_getObjectiveFunc(x);
	}
}

/*
	2進数データを表現型に変換
*/
double GA::_binary2Phenotype(std::vector<int> binary)
{
	// カウント変数
	size_t tmp;

	double decimal;
	int place,  numRow;

	decimal		= 0.0;
	place	= 0;
	for (tmp = binary.size(); tmp >= 0; tmp--)
	{
		x = 0.0;
		place = 0.0;
		if (binary[tmp] > 0)
			decimal += pow(2.0, (double)place);
		place += 1;
	}

	return decimal * 1.0 / (pow(2.0, (double)this->geneLength) - 1.0);
}

double GA::_getObjectiveFunc(double x)
{
	return sin(3.0*x) + 0.5*sin(9.0*x) + sin(15.0*x + 50.0);
}

void GA::output(int generation)
{
	size_t tmp_column;

	std::cout << "第" << generation << "世代" << std::endl;
	for (tmp_column = 0; tmp_column < this->geneLength; tmp_column++)
	{
		x = this->_binary2Phenotype(allIndividual[tmp_column]);
		std::cout << tmp_column << "\t" << x << "\t" << this->_getObjectiveFunc(x) << std::endl;
	}
}

void GA::uniformCrossover()
{
	// カウント変数
	size_t tmp;

	int parent1, parent2;
	std::vector<int> maskPattern, parent1Individual, parent2Individual, child1Individual, child2Individual;
	
	// 親を選択
	do
	{
		parent1 = this->selectIndividual();
		parent2 = this->selectIndividual();
	} while (parent1 == parent2);
	parent1Individual.push_back(this->allIndividual[parent1]);
	parent2Individual.push_back(this->allIndividual[parent2]);

	// マスクパターンを生成
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<int> randomValue(0, 1);
	for (tmp = 0; tmp < this->geneLength; tmp++)
	{
		maskPattern.push_back(randomValue(mt));
	}

	// 交叉
	for (tmp = 0; tmp < this->geneLength; tmp++)
	{
		if (maskPattern[tmp] == 0)
		{
			childIndividual.push_back(parent1Individual[tmp]);
			child2Individual.push_back(parent2Individual[tmp]);
		}
		else if (maskPattern[tmp] == 1)
		{
			childIndividual.push_back(parent2Individual[tmp]);
			child2Individual.push_back(parent1Individual[tmp]);
		}
	}

	// 世代交代
	this->allIndividual[parent1].swap(child1Individual)
	this->allIndividual[parent2].swap(child2Individual);
}

/*
	適応度比例戦略から固体を選択し，その番号を返す．
	失敗した場合は-1を返す．
*/
int GA::selectIndividual()
{
	// カウント変数
	size_t tmp_column;

	std::vector<std::vector<int>> shuffledAllIndividual	= this->allIndividual;
	std::vector<double> shuffledFitness, fitnessRatio, selectedIndividual;
	double sumFitness = 0.0, x;

	// 固体集団をシャッフルして順番を変えずに適応度を求める
	// @MEMO:２重配列のときのsize()の挙動
	std::shuffle(shuffledAllIndividual.begin(), shuffledAllIndividual.end(), std::mt19937());
	for (tmp_column = 0; tmp_column < shuffledAllIndividual.size(); tmp_column++)
	{
		x = this->_binary2Phenotype(shuffledAllIndividual[tmp_column]);
		shuffledFitness.push_back(this->_getObjectiveFunc(x));
		sumFitness += this->_getObjectiveFunc(x);
	}

	// 適応度の比
	for (tmp_column = 0; tmp_column < shuffledAllIndividual.size(); tmp_column++)
	{
		fitnessRatio.push_back(shuffledFitness[tmp_column] / sumFitness);
	}

	// 適応度から確率的に固体を選択
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);
	int selectedFlg = 0;
	do
	{
		for (tmp_column = 0; tmp_column < fitnessRatio.size(); tmp_column++)
		{
			if (fitnessRatio[tmp_column] > randomValue(mt))	// 毎回ランダム値が生成される	
			{
				selectedIndividual.push_back(shuffledAllIndividual[tmp_column]);
				selectedFlg = 1;
				break;
			}
		}
	} while (selectedFlg == 0);

	// 選ばれた個体の番号を返す
	int individualNum;
	bool correspond = false;
	for (individualNum = 0; individualNum < this->population; individualNum++)
	{
		correspond = std::equal(shuffledAllIndividual[tmp_column].cbegin(), shuffledAllIndividual[tmp_column].cend(), selectedIndividual.cbegin());
		if (correspond)
		{
			return individualNum;
		}
	}

	return -1;	// 固体が存在しなかった場合
}

/*
	ランキング法をもとに選択淘汰を行う．
*/
void GA::selectRanking()
{
	// カウント変数
	size_t tmp, tmp1, tmp2;

	// ソート順を記憶するためのキー
	std::vector<int> key;
	for (tmp = 0; tmp < this->population; tmp++)
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

				tmpKey = key[tmp1];
				key[tmp1] = key[tmp2];
				key[tmp2] = tmpKey;
			}
		}
	}

	// 個体集団を適応度の降順に並べる
	for (tmp = 0; tmp < this->geneLength; tmp++)
	{
		this->allIndividual[tmp] = this->allIndividual[key[tmp]];
	}

	this->_assignAllIndividual();
}

void GA::_assignAllIndevidual()
{
	// 【次回】割り当て関数から（本を参考に）
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

	for (tmp_row = 0; tmp_row < this->geneLength; tmp_row++)
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