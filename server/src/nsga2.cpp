#include "../include/nsga2.h"
#include "../include/expfit.h"
#include "../include/common.h"
#include "../include/ga_individual.h"


NSGA2::NSGA2(std::size_t pop, std::size_t iter)
:_popSize(pop),
_iterations(iter)
{}

NSGA2::~NSGA2()
{
	std::vector<GAIndividual>().swap(_finalPops);
}

/**
 * @fn NSGA2を実行する
 * @param MomentEq meq モーメント方程式の情報
 */
int NSGA2::run(MomentEq *meq)
{
	std::size_t i, t;

	// 突然変異と交叉のパラメータ
	double crossProb = 0.85;             // crossover probability
	//double flipProb = 1. / Common::NUM_OF_MOMENTEQ;   // mutation probability
	double flipProb = 0.01;   // mutation probability

	// 変数の定義域を設定
	std::vector<double> lower(Common::NUM_OF_PARAM), upper(Common::NUM_OF_PARAM);
	this->_setValueRange(lower, upper);

	// 親子の個体群を定義
	PopulationMOO parents(_popSize, ChromosomeT< double >(Common::NUM_OF_PARAM));
	PopulationMOO offsprings(_popSize, ChromosomeT< double >(Common::NUM_OF_PARAM));

	// 最小化のタスク
	parents.setMinimize();
	offsprings.setMinimize();

	// 最終的なアーカイ集団
	_archive.setMaxArchive(_popSize);
	_archive.minimize();

	// 親個体群の初期化
	for (i = 0; i < parents.size(); ++i)
	   dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize(0, 1.);
   
	// 目的関数の個数をセット
	parents   .setNoOfObj(meq->getObjNum());
	offsprings.setNoOfObj(meq->getObjNum());

    // 初期個体を生成
	std::vector<double> func(meq->getObjNum());
	ChromosomeT< double > dblchrom;
    for (i = 0; i < parents.size(); ++i) {
 		// モーメント方程式を解く
		dblchrom    = dynamic_cast< ChromosomeT< double > &>(parents[ i ][ 0 ]);
		func = meq->expb_f(dblchrom);
		// 目的関数
		parents[ i ].setMOOFitnessValues(func);
	}

	cout << "NSGA2: start (Number of object:" << meq->getObjNum() << ")" << endl;
	for (t = 1; t <= _iterations; ++t) {
		cout << "generation: " << t << endl;

		// 親個体を子個体にコピー
		offsprings.selectBinaryTournamentMOO(parents);

		// recombine by crossing over two parents
		for (i = 0; i < offsprings.size(); i += 2) {
			if (Rng::coinToss(crossProb)) {
				dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]).SBX(dynamic_cast< ChromosomeT< double > &>(offsprings[ i+1 ][ 0 ]), lower, upper, 20., 0.5);
			}
		}

		// flipping bitsによって突然変異
		for (i = 0; i < offsprings.size(); ++i) {
			dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]).mutatePolynomial(lower, upper, 20., flipProb);
		}

		// 個体群の評価
		for (i = 0; i < parents.size(); ++i) {
			// モーメント方程式を解く
			dblchrom    = dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]);
			func = meq->expb_f(dblchrom);
			// 目的関数
			offsprings[ i ].setMOOFitnessValues(func);
		}

		// 選択
		parents.selectCrowdedMuPlusLambda(offsprings);

		// 10世代おきに劣解をアーカイブから除く
		if (!(t % 10)) {
			_archive.cleanArchive();
			for (i = 0; i < parents.size(); ++i) {
				_archive.addArchive(parents[ i ]);
			}
			_archive.nonDominatedSolutions();
		}
	} // 繰り返し
	cout << "NSGA2: end\n" << endl;

	// data output
	_archive.cleanArchive();
	for (i = 0; i < _popSize; ++i) {
		_archive.addArchive(parents[ i ]);
	}

	_archive.nonDominatedSolutions();

	this->_saveArchive(_archive);
	cout    << "size of the archive: "  << _archive.size() << endl << endl;

	_archive.cleanArchive();

	return EXIT_SUCCESS;
}

/**
 * @fn 変数の範囲を指定する
 * @param vector<double> &lower 下限（領域確保済み）
 * @param vector<double> &upper 上限（領域確保済み）
 */
void NSGA2::_setValueRange(std::vector<double> &lower, std::vector<double> &upper)
{
	lower[0]	= 0.;	upper[0]	= 1.;	// a
	lower[1]	= -3.;	upper[1]	= 3.;	// mu1
	lower[2]	= -2.;	upper[2]	= 2.;	// mu2
	lower[3]	= 0.;	upper[3]	= 1.5;	// sigma11
	lower[4]	= 0.;	upper[4]	= 1.;	// sigma12
	lower[5]	= 0.;	upper[5]	= 1.5;	// sigma21
	lower[6]	= 0.;	upper[6]	= 1.;	// sigma22
	lower[7]	= -1.;	upper[7]	= 1.;	// kappa1
	lower[8]	= -1.;	upper[8]	= 1.;	// kappa2
	lower[9]	= -1.;	upper[9]	= 1.;	// kappa3
}

/**
 * @fn 解の個体群を取得
 * @return vector<GAIndividual> _finalPops 最終的な個体群
 */
std::vector<GAIndividual> NSGA2::getFinalPops()
{
	return _finalPops;
}

/**
 * @fn アーカイブの情報を_finalPopsへ格納
 * @param ArchiveMOO &archive アーカイブ集団
 */
void NSGA2::_saveArchive(ArchiveMOO &archive)
{
	std::size_t i, ii;
	std::size_t no = archive.size();
	std::size_t noOfObj;
	if (no > 0)
		noOfObj = archive.readArchive(0).getNoOfObj();
	else
		noOfObj = 0;

	IndividualMOO individual;
	ChromosomeT< double > chrom;

	_finalPops.resize(no);
	for (i = 0; i < no; ++i)
	{
		individual.operator=(archive.readArchive(i));
		chrom   = dynamic_cast< ChromosomeT< double > &>(individual[0]);
		// モーメント値
		std::vector<double> m;
		m.resize(Common::NUM_OF_MOMENT);
		MomentEq::getMomentFromParameter(chrom, m);
		_finalPops[i].mValue.resize(m.size());
		for (ii = 0; ii < m.size(); ++ii) {
			_finalPops[i].mValue[ii]	= m[ii];
		}
		// 目的関数値
		_finalPops[i].oValue.resize(noOfObj);
		for (ii = 0; ii < noOfObj; ++ii) {
			_finalPops[i].oValue[ii] = archive.readArchive(i).getMOOFitness(ii);
		}
		// パラメータ値
		_finalPops[i].pValue.resize(chrom.size());
		for (ii = 0; ii < chrom.size(); ++ii) {
			_finalPops[i].pValue[ii] = chrom[ii];
		}
	}
}
