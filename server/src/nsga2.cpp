#include "../include/nsga2.h"
#include "../include/expfit.h"
#include "../include/common.h"
#include "../include/ga_individual.h"


NSGA2::NSGA2(int pop, int iter)
{
	_popSize      = pop;  // 120
	_iterations   = iter;  // 500
}

NSGA2::~NSGA2()
{
	std::vector<GAIndividual>().swap(_finalPops);
}

/**
 * @fn NSGA2を実行する
 * @param vecto<double> &dF 入力や系のパラメータ
 */
int NSGA2::run(double *dF)
{
	unsigned int i, ii, t;

	// 突然変異と交叉のパラメータ
	double crossProb = 0.85;             // crossover probability
	double flipProb = 1. / Common::NUM_OF_MOMENTEQ;   // mutation probability
	// double flipProb = 0.01;   // mutation probability

	// 変数の定義域を設定
	std::vector<double> lower(Common::NUM_OF_PARAM), upper(Common::NUM_OF_PARAM);
	this->_setValueRange(lower, upper);

	// 親子の個体群を定義
	PopulationMOO parents(_popSize, ChromosomeT< double >(Common::NUM_OF_PARAM));
	PopulationMOO offsprings(_popSize, ChromosomeT< double >(Common::NUM_OF_PARAM));

	// 最小化のタスク
	parents   .setMinimize();
	offsprings.setMinimize();

	// 最終的なアーカイブ集団
	_archive.setMaxArchive(_popSize);
	_archive.minimize();

	// 親個体群の初期化
	for (i = 0; i < parents.size(); ++i)
	   dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize(0, 1.);
	
	ChromosomeT< double > dblchrom;

	gsl_vector *x	= gsl_vector_alloc(Common::NUM_OF_PARAM);
	gsl_vector *f	= gsl_vector_alloc(Common::NUM_OF_MOMENTEQ);
	gsl_matrix *popData	= gsl_matrix_alloc(Common::NUM_OF_MOMENTEQ, _popSize);

	// 低次元化用の個体群を生成
	for (i = 0; i < _popSize; ++i) {
		// gsl_vectorに変換
		dblchrom    = dynamic_cast< ChromosomeT< double > &>(parents[ i ][ 0 ]);
		for (ii = 0; ii < Common::NUM_OF_PARAM; ++ii) {
			gsl_vector_set(x, ii, dblchrom[ ii ]);
		}
		
		// モーメント方程式を解いてgsl_vectorに変換
		MomentEq::expb_f(x, dF, f);
		for (ii = 0; ii < Common::NUM_OF_MOMENTEQ; ++ii) {
			gsl_matrix_set(popData, ii, i, gsl_vector_get(f, ii));
		}
	}

	// 相関行列（RR^T）
	gsl_matrix *corMtrx;
	corMtrx	= this->_correlationMatrix(popData);
	gsl_matrix_free(popData);
	gsl_matrix *tCorMtrx	= gsl_matrix_alloc(corMtrx->size1, corMtrx->size2);
	gsl_matrix_transpose_memcpy(tCorMtrx, corMtrx);
	gsl_matrix *cor2Mtrx	= gsl_matrix_alloc(corMtrx->size1, corMtrx->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, corMtrx, tCorMtrx, 0.0, cor2Mtrx);
	gsl_matrix_free(tCorMtrx);

	// 固有値解析
	gsl_vector *eval	= gsl_vector_alloc(cor2Mtrx->size1);
	gsl_matrix *evec	= gsl_matrix_alloc(cor2Mtrx->size1, cor2Mtrx->size2);
	gsl_eigen_symmv_workspace *w	= gsl_eigen_symmv_alloc(cor2Mtrx->size1);
	gsl_eigen_symmv(corMtrx, eval, evec, w);	// RR^Tで固有値解析
	gsl_matrix_free(corMtrx);
	gsl_matrix_free(cor2Mtrx);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);	// 寄与率が降順になるようソート

	// 低次元化
	// std::vector<unsigned int> selectedObj = this->_dimensionReduction(eval, evec);
	std::vector<unsigned int> selectedObj{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	// 目的関数の個数をセット
	parents   .setNoOfObj(selectedObj.size());
	offsprings.setNoOfObj(selectedObj.size());

	// 初期個体を生成
	std::vector<double> func(selectedObj.size());
	for (i = 0; i < parents.size(); ++i) {
		for (ii = 0; ii < selectedObj.size(); ++ii) {
			func[ ii ]    = gsl_vector_get(f, selectedObj[ii]);
		}
		// 目的関数
		parents[ i ].setMOOFitnessValues(func);
	}

	cout << "NSGA2: start" << endl;
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
			// gsl_vectorに変換
			dblchrom    = dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]);
			for (ii = 0; ii < Common::NUM_OF_PARAM; ++ii) {
				gsl_vector_set(x, ii, dblchrom[ ii ]);
			}

			// モーメント方程式を解いてgsl_vectorに変換
			MomentEq::expb_f(x, dF, f);
			for (ii = 0; ii < selectedObj.size(); ++ii) {
				func[ ii ]    = gsl_vector_get(f, selectedObj[ii]);
			}

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
	gsl_vector_free(x);
	gsl_vector_free(f);

	return EXIT_FAILURE;
}

/**
 * M(目的関数の数)×N(個体の数)行列から相関行列(M×M)を求める
 * @param gsl_matrix *m M(目的関数の数)×N(個体の数)行列
 * @return gsl_matrix *cor 相関行列
 */
gsl_matrix *
NSGA2::_correlationMatrix(const gsl_matrix *m)
{
	unsigned int i, ii, iii;
	double covValue, corValue;
	gsl_matrix *cov	= gsl_matrix_alloc(m->size1, m->size1);
	gsl_matrix *cor	= gsl_matrix_alloc(m->size1, m->size1);

	// 分散共分散行列の作成
	for (i = 0; i < m->size1; ++i) {
		// i番目の列
		gsl_vector_view rVec1	= gsl_matrix_row(const_cast<gsl_matrix*>(m), i);
		for (ii = 0; ii < m->size1; ++ii) {
			// ii番目の列
			gsl_vector_view rVec2	= gsl_matrix_row(const_cast<gsl_matrix*>(m), ii);
			covValue	= 0.;
			for (iii = 0; iii < m->size2; ++iii) {
				// 分散共分散行列の作成
				covValue	+= gsl_vector_get(&rVec1.vector, iii) * gsl_vector_get(&rVec2.vector, iii);
			}
			gsl_matrix_set(cov, i, ii, covValue/(double)(m->size1 - 1));
		}
	}

	// 相関行列の作成
	for (i = 0; i < m->size1; ++i) {
		for (ii = 0; ii < m->size1; ++ii) {
			if (gsl_matrix_get(cov, i, i) == 0 || gsl_matrix_get(cov, ii, ii) == 0) {
				corValue	= 0.;
			} else {
				corValue	= gsl_matrix_get(cov, i, ii) / sqrt(gsl_matrix_get(cov, i, i)*gsl_matrix_get(cov, ii, ii));
			}
			gsl_matrix_set(cor, i, ii, corValue);
		}
	}

	gsl_matrix_free(cov);
	return cor;
}

/**
 * 目的関数を低次元化する
 * @param gsl_vecor *eval 固有値
 * @param gsl_matrix *evec 固有ベクトル
 * @return vector<unsigned int> rObj 低次元化された目的関数
 */
std::vector<unsigned int>
NSGA2::_dimensionReduction(gsl_vector *eval, gsl_matrix *evec)
{
	unsigned int i;

	// 寄与率を計算
	double sumRate	= 0.;
	for (i = 0; i < eval->size; ++i) {
		sumRate	+= gsl_vector_get(eval, i);
	}	
	gsl_vector *contRate	= gsl_vector_alloc(eval->size);
	for (i = 0; i < eval->size; ++i) {
		gsl_vector_set(contRate, i, gsl_vector_get(eval, i) / sumRate);
	}

	unsigned int cnt = 0;	// 調べる主成分番号
	double cumlCont	= 0.;	// 累積寄与率
	std::vector<unsigned int> rObj, old_rObj;	// 低次元化した目的関数，旧低次元化した目的関数
	rObj.reserve(eval->size*3);

	double mostPos = 0., mostNeg = 0., mostAbs = 0.;	// （固有ベクトル中の）最大の正数，最小の負数，最大の絶対値
	int mostPosNum = -1, mostNegNum = -1, mostAbsNum = -1;	// 最大正数，最小負数，最大絶対値の目的関数
	
	do{
		if (cnt == 0) {
			mostPos = 0.; mostNeg = 0.;
			mostPosNum = -1; mostNegNum = -1;
			// 第一主成分に対して，寄与が一番小さな負数と一番大きな正数の目的関数を選ぶ
			for (i = 0; i < eval->size; ++i) {
				if (mostPos < gsl_matrix_get(evec, cnt, i))	{
					mostPos	= gsl_matrix_get(evec, cnt, i);
					mostPosNum	= i;
				}
				if (mostNeg > gsl_matrix_get(evec, cnt, i))	{
					mostNeg	= gsl_matrix_get(evec, cnt, i);
					mostNegNum	= i;
				}
			}
			if (mostPosNum == -1 || mostNegNum == -1) {
				std::cout << "ERROR: Not exist most positive or most negative objectives in the first principal component." << std::endl;
				exit(1);
			}
			rObj.push_back(mostPosNum);
			rObj.push_back(mostNegNum);
		} else {
			// if (gsl_vector_get(contRate, cnt) <= 0.1) {
			if (gsl_vector_get(eval, cnt) <= 0.1) {
				mostAbs = 0.;
				mostAbsNum = -1;
				for (i = 0; i < eval->size; ++i) {
					if (mostAbs < gsl_matrix_get(evec, cnt, i)) {
						mostAbs = gsl_matrix_get(evec, cnt, i);
						mostAbsNum = i;
					}
				}
				rObj.push_back(mostAbsNum);
			} else if(cumlCont <= 0.95) {
				gsl_vector_const_view rowInEvec = gsl_matrix_const_row(evec, cnt);
				if (gsl_vector_ispos(&rowInEvec.vector)) {
					mostPosNum	= gsl_vector_max_index(&rowInEvec.vector);
					rObj.push_back(mostPosNum);
				} else if (gsl_vector_isneg(&rowInEvec.vector)) {
					rObj.resize(Common::NUM_OF_MOMENTEQ);
					for (i = 0; i < eval->size; ++i) {
						rObj[i] = i;
					}
				} else {
					// 固有ベクトルは正負を含んでいる
					mostPos = gsl_vector_max(&rowInEvec.vector);
					mostPosNum = gsl_vector_max_index(&rowInEvec.vector);
					mostNeg = gsl_vector_min(&rowInEvec.vector);
					mostNegNum = gsl_vector_min_index(&rowInEvec.vector);
					if (mostPos < std::abs(mostNeg)) {
						if (mostPos >= 0.9*std::abs(mostNeg)) {
							rObj.push_back(mostPosNum);
							rObj.push_back(mostNegNum);
						} else {
							rObj.push_back(mostNegNum);
						}
					} else {
						if (mostPos >= 0.8*std::abs(mostNeg)) {
							rObj.push_back(mostPosNum);
							rObj.push_back(mostNegNum);
						} else {
							rObj.push_back(mostPosNum);
						}
					}
				}
			}
		}
		// rObjの重複削除
		std::sort(rObj.begin(), rObj.end());
		rObj.erase( std::unique(rObj.begin(), rObj.end()), rObj.end() );

		cumlCont	+= gsl_vector_get(contRate, cnt);
		if (cnt != 0) {
			if (old_rObj == rObj) break;
		}
		old_rObj = rObj;
		cnt += 1;
	} while(1);

	gsl_vector_free(contRate);
	return rObj;
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
	unsigned int i, ii;
	unsigned int no = archive.size();
	unsigned int noOfObj;
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