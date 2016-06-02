#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/TestFunction.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

using namespace std;

void saveArchiveInFile(char *filename, ArchiveMOO &archive)
{
    /********** 後で消す **********/
    // 境界値
    double   min = 0;
    double   max = 1;
    const Interval rangeOfValues(max, min);
    unsigned numOfBits   = 20;
    bool useGrayCode     = true;
    /********** 後で消す **********/

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
        chrom.decodeBinary(individual[0], rangeOfValues, numOfBits, useGrayCode);
        for (unsigned k = 0; k < chrom.size(); ++k)
        {
            ofs << chrom[k] << " " << std::flush;
        }

        ofs << std::endl;
    }
    ofs.close();
}

// main program
int main(int argc, char* argv[])
{
    unsigned i, ii, k, t;

    unsigned seed        = 0;
    unsigned dimension   = 10;  // problem dimension
    unsigned numOfBits   = 20; 
    unsigned popSize     = 120; // population size
    unsigned iterations  = 500; // maximum number of generations
    bool useGrayCode     = true;

    // 境界値
    double   min = 0;
    double   max = 1;

    // 突然変異と交叉のパラメータ
    double   crossProb = 0.9;            // crossover probability
    double   flipProb = 1. / dimension; // mutation probability

    // 変数の定義域を設定
    // TODO:変数ごとに定義域が異なるものへ対応
    const Interval rangeOfValues(max, min);
    double f1 = 0., f2 = 0.;

    // 親個体と子個体
    double *PF1 = new double[(unsigned)popSize];
    double *PF2 = new double[(unsigned)popSize];
    double *OF1 = new double[(unsigned)popSize];
    double *OF2 = new double[(unsigned)popSize];

    // ランダム値生成器
    Rng::seed(seed);

    // 親子の個体群を定義
    PopulationMOO parents(popSize, ChromosomeT< bool >(dimension*numOfBits));
    PopulationMOO offsprings(popSize, ChromosomeT< bool >(dimension*numOfBits));

    // 最小化のタスク
    parents   .setMinimize();
    offsprings.setMinimize();

    // 目的関数の個数をセット
    parents   .setNoOfObj(2);
    offsprings.setNoOfObj(2);

    // 遺伝子型を表現型にデコードした時の一時保存用
    ChromosomeT< double > dblchrom;

    // 親個体群の初期化
    for (i = 0; i < parents.size(); ++i)
       dynamic_cast< ChromosomeT< bool >& >( parents[ i ][ 0 ] ).initialize();
        // ((ChromosomeT< double >*)&parents[ i ][ 0 ])->initialize(minInit, maxInit);
    
    // 個体群の評価
    for (i = 0; i < parents.size(); ++i) {
        dblchrom.decodeBinary(parents[ i ][ 0 ], rangeOfValues, numOfBits, useGrayCode);

        // 目的関数
        f1 = ZDT4FF1(dblchrom);
        f2 = ZDT4FF2(dblchrom);
        
        parents[ i ].setMOOFitnessValues(f1, f2);
        PF1[ i ] = f1;
        PF2[ i ] = f2;
    }

    // iterate
    for (t = 1; t <= iterations; ++t) {
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
            dblchrom.decodeBinary(offsprings[ i ][ 0 ], rangeOfValues, numOfBits, useGrayCode);

            // 目的関数
            f1 = ZDT4FF1(dblchrom);
            f2 = ZDT4FF2(dblchrom);
        
            offsprings[ i ].setMOOFitnessValues(f1, f2);
            OF1[ i ] = f1;
            OF2[ i ] = f2;
        }

        // 選択
        parents.selectCrowdedMuPlusLambda(offsprings);
        for (k = 0; k < parents.size(); k++) {
            PF1[ k ] = parents[ k ].getMOOFitness(0);
            PF2[ k ] = parents[ k ].getMOOFitness(1);
        }
    } // 繰り返し

    // data output
    ArchiveMOO archive(popSize);
    archive.minimize();
    for (ii = 0; ii < popSize; ++ii) {
        archive.addArchive(parents[ ii ]);
    }
    archive.nonDominatedSolutions();

    char filename[10000];
    sprintf(filename, "nsga2example-%d.out", 1);
    saveArchiveInFile(filename, archive);

    cout    << "size of the archive: "  << archive.size()
    << ", filename of archive: " << filename << endl << endl;

    delete [] PF1;
    delete [] PF2;
    delete [] OF1;
    delete [] OF2;

    return EXIT_SUCCESS;
}
