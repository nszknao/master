#include "../include/nsga2.h"
#include "../include/expfit.h"

// コンストラクタ
NSGA2::NSGA2(int pop, int bits, bool gray, int iter) {
    // 境界値
    _max  = 0;
    _min  = 1;
    _popSize      = pop;  // 120
    _numOfBits    = bits; // 20
    _useGrayCode  = gray; // true
    _iterations   = iter;  // 500
}

/*
    NSGA2を実行する
*/
int NSGA2::run(nsga2_function * function)
{
    unsigned i, ii, k, t, ret;

    unsigned seed   = 0;

    // 突然変異と交叉のパラメータ
    double crossProb = 0.9;            // crossover probability
    double flipProb = 1. / function->n; // mutation probability

    // 変数の定義域を設定
    // TODO:各パラメータで設定できるように
    const Interval rangeOfValues(_max, _min);

    // 親個体と子個体
    std::vector< std::vector<double> > PF(function->n, std::vector<double>(_popSize));
    std::vector< std::vector<double> > OF(function->n, std::vector<double>(_popSize));

    // ランダム値生成器
    Rng::seed(seed);

    // 親子の個体群を定義
    PopulationMOO parents(_popSize, ChromosomeT< bool >(function->n*_numOfBits));
    PopulationMOO offsprings(_popSize, ChromosomeT< bool >(function->n*_numOfBits));

    // 最小化のタスク
    parents   .setMinimize();
    offsprings.setMinimize();

    // 目的関数の個数をセット
    parents   .setNoOfObj(function->n);
    offsprings.setNoOfObj(function->n);

    // 遺伝子型を表現型にデコードした時の一時保存用
    ChromosomeT< double > dblchrom;

    // 親個体群の初期化
    // TODO:初期値は自分で指定
    for (i = 0; i < parents.size(); ++i)
       dynamic_cast< ChromosomeT< bool >& >( parents[ i ][ 0 ] ).initialize();

    // 目的関数値を保存する
    std::vector<double> func(function->n);
    
    gsl_vector *x, *f;
    // 個体群の評価
    for (i = 0; i < parents.size(); ++i) {
        dblchrom.decodeBinary(parents[ i ][ 0 ], rangeOfValues, _numOfBits, _useGrayCode);

        // 変数をgsl_vectorにセット
        for (ii = 0; ii < function->p; ++ii)
            gsl_vector_set(x, t, dblchrom[ ii ]);
        // 方程式
        ret = MomentEq::expb_f(x, function->params, f);
        for (ii = 0; ii < function->n; ++i)
            func[ ii ]    = gsl_vector_get(f, ii);

        // 目的関数
        parents[ i ].setMOOFitnessValues(func);

        // 親個体
        for (ii = 0; ii < function->n; ++ii)
            PF[ i ][ ii ] = func[ ii ];
    }

    // iterate
    for (t = 1; t <= _iterations; ++t) {
        cout << "Generation: " << t << endl;

        // 親個体を子個体にコピー
        offsprings.selectBinaryTournamentMOO(parents);

        // recombine by crossing over two parents
        for (i = 0; i < offsprings.size(); i += 2)
            if (Rng::coinToss(crossProb))
                offsprings[ i ][ 0 ].crossoverUniform(offsprings[ i+1 ][ 0 ]);

        // flipping bitsによって突然変異
        for (i = 0; i < offsprings.size(); ++i)
            dynamic_cast< ChromosomeT< bool >& >(offsprings[ i ][ 0 ]).flip(flipProb);

        // 個体群の評価
        for (i = 0; i < parents.size(); ++i) {
            dblchrom.decodeBinary(offsprings[ i ][ 0 ], rangeOfValues, _numOfBits, _useGrayCode);

            // 変数をgsl_vectorにセット
            for (ii = 0; ii < function->p; ++ii)
                gsl_vector_set(x, t, dblchrom[ ii ]);
            // 方程式
            ret = MomentEq::expb_f(x, function->params, f);
            for (ii = 0; ii < function->n; ++ii)
                func[ ii ]    = gsl_vector_get(f, ii);

            // 目的関数
            offsprings[ i ].setMOOFitnessValues(func);

            // 子個体
            for (ii = 0; ii < function->n; ++ii)
                OF[ i ][ ii ] = func[ ii ];
        }

        // 選択
        parents.selectCrowdedMuPlusLambda(offsprings);
        for (k = 0; k < parents.size(); k++) {
            for (ii = 0; ii < function->n; ++ii) {
                PF[ k ][ ii ]   = parents[ k ].getMOOFitness(ii);
            }
        }
    } // 繰り返し

    // data output
    ArchiveMOO archive(_popSize);
    archive.minimize();
    for (ii = 0; ii < _popSize; ++ii) {
        archive.addArchive(parents[ ii ]);
    }
    archive.nonDominatedSolutions();

    char filename[10000];
    sprintf(filename, "nsga2example-%d.out", 1);
    this->_saveArchiveInFile(filename, archive);

    cout    << "size of the archive: "  << archive.size()
    << ", filename of archive: " << filename << endl << endl;

    return EXIT_SUCCESS;
}

/*
    アーカイブの情報をファイルに書き込む．
    目的関数値1 目的関数値2 ... | パラメータ値1 パラメータ値2 ...
    @param *filename ファイルポインタ
    @param &archive アーカイブ集団
*/
void NSGA2::_saveArchiveInFile(char *filename, ArchiveMOO &archive)
{
    // TODO:各パラメータで設定できるように
    const Interval rangeOfValues(_max, _min);

    unsigned no = archive.size();
    unsigned noOfObj;
    if (no > 0)
        noOfObj = archive.readArchive(0).getNoOfObj();
    else
        noOfObj = 0;

    IndividualMOO individual;
    ChromosomeT< double > chrom;

    double f;
    ofstream ofs(filename);
    for (unsigned i = 0; i < no; ++i)
    {
        // 目的関数値
        for (unsigned j = 0; j < noOfObj; ++j)
        {
            f   = archive.readArchive(i).getMOOFitness(j);
            ofs << f << " " << std::flush;
        }

        ofs << "| ";

        // パラメータ値
        individual.operator=(archive.readArchive(i));
        chrom.decodeBinary(individual[0], rangeOfValues, _numOfBits, _useGrayCode);
        for (unsigned k = 0; k < chrom.size(); ++k)
        {
            ofs << chrom[k] << " " << std::flush;
        }

        ofs << std::endl;
    }
    ofs.close();
}