#include "../include/nsga2.h"
#include "../include/expfit.h"

NSGA2::NSGA2(int pop, int iter) {
	_popSize      = pop;  // 120
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
	std::vector<double> lower(params->p), upper(params->p);
	this->_setValueRange(lower, upper);

	// 親個体と子個体
	std::vector< std::vector<double> > PF(params->n, std::vector<double>(_popSize));
	std::vector< std::vector<double> > OF(params->n, std::vector<double>(_popSize));

	// ランダム値生成器
	Rng::seed(seed);

	// 親子の個体群を定義
	PopulationMOO parents(_popSize, ChromosomeT< double >(params->p));
	PopulationMOO offsprings(_popSize, ChromosomeT< double >(params->p));

	// 最小化のタスク
	parents   .setMinimize();
	offsprings.setMinimize();

	// 目的関数の個数をセット
	parents   .setNoOfObj(params->n);
	offsprings.setNoOfObj(params->n);

	// 最終的なアーカイブ集団
	_archive.setMaxArchive(_popSize);
	_archive.minimize();


	// 親個体群の初期化
	for (i = 0; i < parents.size(); ++i)
	   dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize(0, 1.);
	
	ChromosomeT< double > dblchrom;

	// 目的関数値を保存する
	std::vector<double> func(params->n);

	gsl_vector *x, *f;
	x   = gsl_vector_alloc(params->p);
	f   = gsl_vector_alloc(params->n);
	// 個体群の評価
	for (i = 0; i < parents.size(); ++i) {
		// gsl_vectorに変換
		dblchrom    = dynamic_cast< ChromosomeT< double > &>(parents[ i ][ 0 ]);
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
		cout << "generation: " << t << endl;

		// 親個体を子個体にコピー
		offsprings.selectBinaryTournamentMOO(parents);

		// recombine by crossing over two parents
		for (i = 0; i < offsprings.size(); i += 2) {
			if (Rng::coinToss(crossProb))
				dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]).SBX(dynamic_cast< ChromosomeT< double > &>(offsprings[ i+1 ][ 0 ]), lower, upper, 20., 0.5);
		}

		// flipping bitsによって突然変異
		for (i = 0; i < offsprings.size(); ++i)
			dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]).mutatePolynomial(lower, upper, 20., flipProb);

		// 個体群の評価
		for (i = 0; i < parents.size(); ++i) {
			// gsl_vectorに変換
			dblchrom    = dynamic_cast< ChromosomeT< double > &>(offsprings[ i ][ 0 ]);
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

		// 10世代おきに劣解をアーカイブから除く
		if (!(t % 10)) {
			_archive.cleanArchive();
			for (i = 0; i < parents.size(); ++i)
				_archive.addArchive(parents[ i ]);
			_archive.nonDominatedSolutions();
		}
	} // 繰り返し
	cout << "NSGA2: end\n" << endl;

	// data output
	_archive.cleanArchive();
	for (i = 0; i < _popSize; ++i)
		_archive.addArchive(parents[ i ]);

	_archive.nonDominatedSolutions();

	cout << _archive.size() << endl;
	this->_saveArchive(_archive);

	cout    << "size of the archive: "  << _archive.size() << endl << endl;

	gsl_vector_free(x);
	gsl_vector_free(f);

	return EXIT_SUCCESS;
}

/**
 * @fn アーカイブの情報をファイルに書き込む．
 * 目的関数値1 目的関数値2 ... | パラメータ値1 パラメータ値2 ...
 * @param string name ファイル名
 * @param ArchiveMOO &archive アーカイブ集団
 */
void NSGA2::saveArchiveInFile(const std::string name)
{
	unsigned int i, ii;
	unsigned int no = _archive.size();
	unsigned int noOfObj;
	if (no > 0)
		noOfObj = _archive.readArchive(0).getNoOfObj();
	else
		noOfObj = 0;

	IndividualMOO individual;
	ChromosomeT< double > chrom;

	double f;
	ofstream ofs(name);
	for (i = 0; i < no; ++i)
	{
		individual.operator=(_archive.readArchive(i));
		chrom   = dynamic_cast< ChromosomeT< double > &>(individual[0]);
		// モーメント値
		std::vector<double> m;
		MomentEq::getMomentFromParameter(chrom, m);
		for (ii = 0; ii < m.size(); ++ii) {
			ofs <<  m[ii] << " " << std::flush;
		}
		// 目的関数値（絶対値）
		for (ii = 0; ii < noOfObj; ++ii) {
			f   = _archive.readArchive(i).getMOOFitness(ii);
			ofs << fabs(f) << " " << std::flush;
		}
		// パラメータ値
		for (ii = 0; ii < chrom.size(); ++ii) {
			ofs << chrom[ii] << " " << std::flush;
		}

		ofs << std::endl;
	}
}

/**
 * @fn 変数の範囲を指定する
 * @param std::vector<double> &lower 下限（領域確保済み）
 * @param std::vector<double> &upper 上限（領域確保済み）
 */
void NSGA2::_setValueRange(std::vector<double> &lower, std::vector<double> &upper)
{
	lower[0]	= 0.;	upper[0]	= 1.;	// a
	lower[1]	= -3.;	upper[1]	= 3.;	// mu1
	lower[2]	= -2.;	upper[2]	= 2.;	// mu2
	lower[3]	= 0.;	upper[3]	= 1.0;	// sigma11
	lower[4]	= 0.;	upper[4]	= 1.0;	// sigma12
	lower[5]	= 0.;	upper[5]	= 1.5;	// sigma21
	lower[6]	= 0.;	upper[6]	= 1.0;	// sigma22
	lower[7]	= -1.;	upper[7]	= 1.;	// kappa1
	lower[8]	= -1.;	upper[8]	= 1.;	// kappa2
	lower[9]	= -1.;	upper[9]	= 1.;	// kappa3
}

std::vector< std::vector<double> > NSGA2::getObjValue()
{
	return _obj;
}

std::vector< std::vector<double> > NSGA2::getPrmValue()
{
	return _prm;
}

/**
 * @fn アーカイブの情報を格納
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
	Common::resize2DemensionalVector(_obj, no, noOfObj);

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
		chrom   = dynamic_cast< ChromosomeT< double > &>(individual[0]);
		if (i == 0)
			Common::resize2DemensionalVector(_prm, no, chrom.size());

		for (ii = 0; ii < chrom.size(); ++ii)
			_prm[i][ii] = chrom[ii];
	}
}
