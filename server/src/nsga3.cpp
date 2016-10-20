#include "../include/nsga3.h"

#include "nsga3cpp/src/problem_base.h"
#include "nsga3cpp/src/alg_nsgaiii.h"
#include "nsga3cpp/src/exp_experiment.h"
#include "nsga3cpp/src/alg_population.h"

#include "nsga3cpp/src/log.h"
#include "nsga3cpp/src/aux_misc.h"

#include "../include/expfit.h"

using namespace std;

NSGA3::~NSGA3()
{
	std::vector<GAIndividual>().swap(_finalPops);
}

/**
 * @fn NSGA3を実行する
 */
int NSGA3::run()
{
	// 設定ファイルの存在確認
	ifstream exp_list("/usr/local/src/master/src/nsga3cpp/explist.ini");
	if (!exp_list) { cout << "We need the explist.ini file." << endl; return 1; }

	string exp_name;
	while (exp_list >> exp_name)
	{
		// 設定ファイル中の実験リストを読み込む
		ifstream exp_ini("/usr/local/src/master/src/nsga3cpp/Experiments/" + exp_name);
		if (!exp_ini) { cout << exp_name << " file does not exist." << endl; continue; }

		// ----- Setup the expriment ------
		CNSGAIII nsgaiii;	// alg_nsgaiii.h
		BProblem *problem = 0;	// problem_base.h

		SetupExperiment(nsgaiii, &problem, exp_ini);	// exp_experiment.h
		// IGD(Inverted Generational Distance)
		// IGDは得られたパレートフロントが，どの程度パレート最適フロントに類似しているかを表す．
		// @link http://www.is.doshisha.ac.jp/academic/papers/pdf/09/2009mthesis/2009luy.pdf
		ofstream IGD_results(nsgaiii.name() + "-" + problem->name() + "-IGD.txt"); // output file for IGD values per run


		// ----- Run the algorithm to solve the designated function -----
		const size_t NumRuns = 1; // 20 is the setting in NSGA-III paper
		for (size_t r=0; r<NumRuns; r+=1)
		{
			srand(r); cout << "Solving " << problem->name() << " ... Run: " << r << endl;

			// --- Solve
			CPopulation solutions;	// alg_population.h
			nsgaiii.Solve(&solutions, *problem);

			// --- Output the result
			string logfname = "/usr/local/src/master/src/nsga3cpp/Results/" + nsgaiii.name() + "-" + problem->name() + "-Run" + IntToStr(r) + ".txt"; // e.g. NSGAIII-DTLZ1(3)-Run0.txt
			SaveScatterData(logfname, solutions, ldAll);	// log.h

			// --- Store the result
			if (r == NumRuns-1) {
				this->_savePopulation(solutions);
			}
		}
		delete problem;

	}// while - there are more experiments to carry out

	return 0;
}

/**
 * @fn 解の個体群を取得
 * @return vector<GAIndividual> _finalPops 最終的な個体群
 */
std::vector<GAIndividual> NSGA3::getFinalPops()
{
	return _finalPops;
}

/**
 * @fn アーカイブの情報を_finalPopsへ格納
 * @param CPopulation &population アーカイブ集団
 */
void NSGA3::_savePopulation(const CPopulation &pop)
{
	unsigned int i, ii;
	
	_finalPops.resize(pop.size());
	for (i = 0; i < pop.size(); ++i) {
		// モーメント値
		std::vector<double> m;
		m.resize(21);
		MomentEq::getMomentFromParameter(pop[i].vars(), m);
		_finalPops[i].mValue.resize(m.size());
		for (ii = 0; ii < m.size(); ++ii) {
			_finalPops[i].mValue[ii]	= m[ii];
		}
		// 目的関数値
		_finalPops[i].oValue.resize(pop[i].objs().size());
		for (ii = 0; ii < pop[i].objs().size(); ++ii) {
			_finalPops[i].oValue[ii]	= pop[i].objs()[ii];
		}
		// パラメータ値
		_finalPops[i].pValue.resize(pop[i].vars().size());
		for (ii = 0; ii < pop[i].vars().size(); ++ii) {
			_finalPops[i].pValue[ii]	= pop[i].vars()[ii];
		}
	}
}
