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
    double crossProb = 0.95;             // crossover probability
//    double flipProb = 1. / Common::NUM_OF_MOMENTEQ;   // mutation probability
    double flipProb = 0.01;   // mutation probability

    std::vector<double> lower = meq->getLowerObj();
    std::vector<double> upper = meq->getUpperObj();

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
    for (i = 0; i < parents.size(); ++i) {
        dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize(lower, upper);
    }
   
    // 目的関数の個数をセット
    parents   .setNoOfObj(meq->getObjNum());
    offsprings.setNoOfObj(meq->getObjNum());

    // 初期個体を生成
    std::vector<double> func(meq->getObjNum()); // temporary addressエラー回避のために必要
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
//        if (!(t % 10)) {
//            _archive.cleanArchive();
//            for (i = 0; i < parents.size(); ++i) {
//                _archive.addArchive(parents[ i ]);
//            }
//            _archive.nonDominatedSolutions();
//        }
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
        std::vector<double> m = MomentEq::getMomentFromParameter(chrom);
        _finalPops[i].mValue.resize(m.size());
        for (ii = 0; ii < m.size(); ++ii) {
            _finalPops[i].mValue[ii]    = m[ii];
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
