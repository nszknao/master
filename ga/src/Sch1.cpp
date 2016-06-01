/*
* <Sch1.cpp>
* Copyright (C) DOLPHIN Project-Team, INRIA Futurs, 2006-2007
* (C) OPAC Team, LIFL, 2002-2007
*
* Abdelhakim Deneche, Arnaud Liefooghe
*
* This software is governed by the CeCILL license under French law and
* abiding by the rules of distribution of free software.  You can  use,
* modify and/ or redistribute the software under the terms of the CeCILL
* license as circulated by CEA, CNRS and INRIA at the following URL
* "http://www.cecill.info".
*
* As a counterpart to the access to the source code and  rights to copy,
* modify and redistribute granted by the license, users are provided only
* with a limited warranty  and the software's author,  the holder of the
* economic rights,  and the successive licensors  have only  limited liability.
*
* In this respect, the user's attention is drawn to the risks associated
* with loading,  using,  modifying and/or developing or reproducing the
* software by the user in light of its specific status of free software,
* that may mean  that it is complicated to manipulate,  and  that  also
* therefore means  that it is reserved for developers  and  experienced
* professionals having in-depth computer knowledge. Users are therefore
* encouraged to load and test the software's suitability as regards their
* requirements in conditions enabling the security of their systems and/or
* data to be ensured and,  more generally, to use and operate it in the
* same conditions as regards security.
* The fact that you are presently reading this means that you have had
* knowledge of the CeCILL license and that you accept its terms.
*
* ParadisEO WebSite : http://paradiseo.gforge.inria.fr
* Contact: paradiseo-help@lists.gforge.inria.fr
*
*/
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <moeo>
#include <es/eoRealInitBounded.h>
#include <es/eoRealOp.h>

#include <math.h>

using namespace std;


/*
    １．目的関数の数を定義
    ２．最小化or最大化を指定
    moeoObjectiveVectorTritsを継承している必要がある．
*/
class Sch1ObjectiveVectorTraits : public moeoObjectiveVectorTraits
{
public:
    static bool minimizing (int i)
    {
        return true;
    }
    static bool maximizing (int i)
    {
        return false;
    }
    static unsigned int nObjectives ()
    {
        return 3;
    }
};


/*
    目的関数ベクトルを定義
    実数値として目的関数ベクトルを用いる．
*/
typedef moeoRealObjectiveVector < Sch1ObjectiveVectorTraits > Sch1ObjectiveVector;


// multi-objective evolving object for the Sch1 problem
/*
    解の定義
    実数値の遺伝子型を変数の数だけ用意する．
*/
class Sch1 : public moeoRealVector < Sch1ObjectiveVector >
{
public:
    Sch1() : moeoRealVector < Sch1ObjectiveVector > (10)
    {}
};


// evaluation of objective functions
/*
    評価する目的関数を定義
    Sch1でテンプレート化されたmoeoEvalFuncを継承している必要がある．
*/
class Sch1Eval : public moeoEvalFunc < Sch1 >
{
public:
    void operator () (Sch1 & _sch1)
    {
        if (_sch1.invalidObjectiveVector())
        {
            Sch1ObjectiveVector objVec;
            double x1 = _sch1[0], x2 = _sch1[1], x3 = _sch1[2], x4 = _sch1[3], x5 = _sch1[4];
            double x6 = _sch1[5], x7 = _sch1[6], x8 = _sch1[7], x9 = _sch1[8], x10 = _sch1[9];
            // 不等式の制約条件
            double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;
            (0 <= x1 && x1 <= 1) ? c1 = 0. : c1 = x1;     (0 <= x2 && x2 <= 1) ? c2 = 0. : c2 = x2;
            (0 <= x3 && x3 <= 1) ? c3 = 0. : c3 = x3;     (0 <= x4 && x4 <= 1) ? c4 = 0. : c4 = x4;
            (0 <= x5 && x5 <= 1) ? c5 = 0. : c5 = x5;     (0 <= x6 && x6 <= 1) ? c6 = 0. : c6 = x6;
            (0 <= x7 && x7 <= 1) ? c7 = 0. : c7 = x7;     (0 <= x8 && x8 <= 1) ? c8 = 0. : c8 = x8;
            (0 <= x9 && x9 <= 1) ? c9 = 0. : c9 = x9;     (0 <= x10 && x10 <= 1) ? c10 = 0. : c10 = x10;
            double g, h, f1, f2, N = 10.;


            f1 = x1;
            g = 1. + 10.*(x2+x3+x4+x5+x6+x7+x8+x9+x10)/((double)N-1.);
            h = 1. - pow(f1/g,0.25) - f1/g*sin(10.*M_PI*f1);
            f2 = g*h;


            objVec[0] = f1;
            objVec[1] = f2;
            objVec[2] = 1000000.*pow(c1,2) + pow(c2,2) + pow(c3,2) + pow(c4,2) + pow(c5,2)
                                + pow(c6,2) + pow(c7,2) + pow(c8,2) + pow(c9,2) + pow(c10,2);
            cout << objVec[0] << " " << objVec[1] << endl;
            // double x = _sch1[0];
            // objVec[0] = pow(x,2);
            // objVec[1] = pow(x-2.,2);
            _sch1.objectiveVector(objVec);
        }
    }
};


int main (int argc, char *argv[])
{
    eoParser parser(argc, argv);  // for user-parameter reading
    eoState state;                // to keep all things allocated

    /* 各種パラメータの定義 */
    // 個体数
    unsigned int POP_SIZE = parser.createParam((unsigned int)(120), "popSize", "Population size",'P',"Param").value();
    // 最大の世代数
    unsigned int MAX_GEN = parser.createParam((unsigned int)(250), "maxGen", "Maximum number of generations",'G',"Param").value();
    // ??（rangeとは）
    double M_EPSILON = parser.createParam(0.01, "mutEpsilon", "epsilon for mutation",'e',"Param").value();
    // ??（交叉確率）
    double P_CROSS = parser.createParam(0.5, "pCross", "Crossover probability",'C',"Param").value();
    // 突然変異率
    double P_MUT = parser.createParam(0.01, "pMut", "Mutation probability",'M',"Param").value();

    // objective functions evaluation
    Sch1Eval eval;

    /* 交叉と突然変異（種類が多く用意されている） */
    eoHypercubeCrossover < Sch1 > xover;
    // eoSegmentCrossover < Sch1 > xover;
    // eoQuadCloneOp < Sch1 > xover;
    eoUniformMutation < Sch1 > mutation (M_EPSILON);

    /* 初期個体の生成 */
    // 境界値の指定（引数は，（変数の数）（最小値）（最大値））
    vector<double> minBounds{0,0,0,0,0,0,0,0,0,0};
    vector<double> maxBounds{1,1,1,1,1,1,1,1,1,1};
    eoRealVectorBounds bounds (minBounds, maxBounds);
    eoRealInitBounded < Sch1 > init (bounds);
    eoPop < Sch1 > pop (POP_SIZE, init);

    // build NSGA-II
    moeoNSGAII < Sch1 > nsgaII (MAX_GEN, eval, xover, P_CROSS, mutation, P_MUT);

    // help ?
    make_help(parser);

    // run the algo
    nsgaII (pop);

    // extract first front of the final population using an moeoArchive (this is the output of nsgaII)
    /*  */
    moeoUnboundedArchive < Sch1 > arch;
    arch(pop);

    // printing of the final archive
    /* 結果を表示 */
    // 表示形式は，（目的関数1の値）（目的関数2の値）・・・（変数の数）（変数1の値）（変数2の値）・・・
    cout << "Final Archive" << endl;
    arch.sortedPrintOn (cout);
    cout << endl;

    return EXIT_SUCCESS;
}
