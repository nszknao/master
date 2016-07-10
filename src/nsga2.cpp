#include "../include/nsga2.h"
#include "../include/expfit.h"

NSGA2::NSGA2(int pop, int bits, bool gray, int iter) {
    // 境界値
    _max  = 0;
    _min  = 1;
    _popSize      = pop;  // 120
    _numOfBits    = bits; // 20
    _useGrayCode  = gray; // true
    _iterations   = iter;  // 500
}

void NSGA2::freeVector()
{
    std::vector< std::vector<double> >().swap(_obj);
    std::vector< std::vector<double> >().swap(_prm);
}

/**
 * @fn NSGA2を実行する
 * @param ParamData* params 入力や系のパラメータ
 */
int NSGA2::run(ParamData* params)
{
    unsigned int i, ii, k, t;

    unsigned seed   = 0;

    // 突然変異と交叉のパラメータ
    double crossProb = 0.9;             // crossover probability
    double flipProb = 1. / params->n;   // mutation probability

    // 変数の定義域を設定
    // TODO:各パラメータで設定できるように
    const Interval rangeOfValues(_max, _min);

    // 親個体と子個体
    std::vector< std::vector<double> > PF(params->n, std::vector<double>(_popSize));
    std::vector< std::vector<double> > OF(params->n, std::vector<double>(_popSize));

    // ランダム値生成器
    Rng::seed(seed);

    // 親子の個体群を定義
    PopulationMOO parents(_popSize, ChromosomeT< bool >(params->p*_numOfBits));
    PopulationMOO offsprings(_popSize, ChromosomeT< bool >(params->p*_numOfBits));

    // 最小化のタスク
    parents   .setMinimize();
    offsprings.setMinimize();

    // 目的関数の個数をセット
    parents   .setNoOfObj(params->n);
    offsprings.setNoOfObj(params->n);

    // 遺伝子型を表現型にデコードした時の一時保存用
    ChromosomeT< double > dblchrom;

    // 親個体群の初期化
    // TODO:初期値は自分で指定
    for (i = 0; i < parents.size(); ++i)
       dynamic_cast< ChromosomeT< bool >& >( parents[ i ][ 0 ] ).initialize();

    // 目的関数値を保存する
    std::vector<double> func(params->n);

    gsl_vector *x, *f;
    x   = gsl_vector_alloc(params->p);
    f   = gsl_vector_alloc(params->n);
    // 個体群の評価
    for (i = 0; i < parents.size(); ++i) {
        // 個体をデコードしてgsl_vectorに変換
        dblchrom.decodeBinary(parents[ i ][ 0 ], rangeOfValues, _numOfBits, _useGrayCode);
        for (ii = 0; ii < params->p; ++ii)
            gsl_vector_set(x, ii, dblchrom[ ii ]);
        
        // モーメント方程式を解いてgsl_vectorに変換
        MomentEq::expb_f(x, params, f);
        for (ii = 0; ii < params->n; ++ii)
            func[ ii ]    = gsl_vector_get(f, ii);

        // 目的関数
        parents[ i ].setMOOFitnessValues(func);

        // 親個体
        for (ii = 0; ii < params->n; ++ii)
            PF[ ii ][ i ] = func[ ii ];
    }

    // iterate
    cout << "NSGA2: start" << endl;
    for (t = 1; t <= _iterations; ++t) {
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
            // 個体をデコードしてgsl_vectorに変換
            dblchrom.decodeBinary(offsprings[ i ][ 0 ], rangeOfValues, _numOfBits, _useGrayCode);
            for (ii = 0; ii < params->p; ++ii)
                gsl_vector_set(x, ii, dblchrom[ ii ]);

            // モーメント方程式を解いてgsl_vectorに変換
            MomentEq::expb_f(x, params, f);
            for (ii = 0; ii < params->n; ++ii)
                func[ ii ]    = gsl_vector_get(f, ii);

            // 目的関数
            offsprings[ i ].setMOOFitnessValues(func);

            // 子個体
            for (ii = 0; ii < params->n; ++ii)
                OF[ ii ][ i ] = func[ ii ];
        }

        // 選択
        parents.selectCrowdedMuPlusLambda(offsprings);
        for (k = 0; k < parents.size(); ++k) {
            for (ii = 0; ii < params->n; ++ii) {
                PF[ ii ][ k ]   = parents[ k ].getMOOFitness(ii);
            }
        }
    } // 繰り返し
    cout << "NSGA2: end\n" << endl;

    // data output
    ArchiveMOO archive(_popSize);
    archive.minimize();
    for (i = 0; i < _popSize; ++i) {
        archive.addArchive(parents[ i ]);
    }
    archive.nonDominatedSolutions();

    this->_saveArchive(archive);
    // this->_saveArchiveInFile((char*)"nsga2example.out", archive);

    cout    << "size of the archive: "  << archive.size() << endl << endl;

    gsl_vector_free(x);
    gsl_vector_free(f);

    return EXIT_SUCCESS;
}

std::vector< std::vector<double> > NSGA2::getObjValue()
{
    return _obj;
}

std::vector< std::vector<double> > NSGA2::getPrmValue()
{
    return _prm;
}

/*
    アーカイブの情報を格納
    @param &archive アーカイブ集団
*/
void NSGA2::_saveArchive(ArchiveMOO &archive)
{
    // TODO:各パラメータで設定できるように
    const Interval rangeOfValues(_max, _min);

    unsigned int i, ii, iii;
    unsigned int no = archive.size();
    unsigned int noOfObj;
    if (no > 0)
        noOfObj = archive.readArchive(0).getNoOfObj();
    else
        noOfObj = 0;
    _obj.resize(no);
    for (i = 0; i < no; ++i) {
        _obj[i].resize(noOfObj);
    }


    IndividualMOO individual;
    ChromosomeT< double > chrom;

    double f;
    for (i = 0; i < no; ++i)
    {
        // 目的関数値
        for (ii = 0; ii < noOfObj; ++ii) {
            f   = archive.readArchive(i).getMOOFitness(ii);
            _obj[i][ii] = f;
        }

        // パラメータ値
        individual.operator=(archive.readArchive(i));
        chrom.decodeBinary(individual[0], rangeOfValues, _numOfBits, _useGrayCode);
        if (i == 0) {
            _prm.resize(no);
            for (iii = 0; iii < no; ++iii) {
                _prm[iii].resize(chrom.size());
            }
        }
        for (ii = 0; ii < chrom.size(); ++ii)
            _prm[i][ii] = chrom[ii];
    }
}

/**
 * @fn ２次元vectorのリサイズ
 */

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