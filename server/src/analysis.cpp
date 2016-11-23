#include "../include/analysis.h"
#include "../include/expfit.h"
#include "../include/nsga2.h"
#include "../include/nsga3.h"
#include "../include/common.h"
#include "../include/ga_individual.h"

const double Analysis::GGD_KAPPA = 2.;	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布

Analysis::Analysis(double arg_l, double arg_b, double arg_a)
:_lambda(arg_l),
_beta2(arg_b),
_alpha(arg_a)
{}

// デストラクタ
Analysis::~Analysis() {}

/**
 * @fn NSGA2でモーメント方程式を解く
 * @param vector<double> &pValue 計算結果のパラメータ値を保存
 * @param vector<double> &oValue 計算結果の目的関数値を保存
 */
int Analysis::GeneticAlgorithm(std::vector<GAIndividual> &pops)
{
	std::cout << "Get analysis solution using Genetic Algorithm.\n" << std::endl;

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / Analysis::GGD_KAPPA)*pow(gsl_sf_gamma(3. / Analysis::GGD_KAPPA), -1.)*this->_beta2);

	// 入力に関するモーメント
	std::vector<double> dF(6);
	dF[0] = 0;
	dF[1] = this->_alpha*Common::S0 + this->_lambda*(1. - this->_alpha)*this->_beta2;
	dF[2] = 0;
	dF[3] = this->_lambda*pow((1. - this->_alpha), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / Analysis::GGD_KAPPA)*pow(gsl_sf_gamma(1. / Analysis::GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = this->_lambda*pow((1. - this->_alpha), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / Analysis::GGD_KAPPA)*pow(gsl_sf_gamma(1. / Analysis::GGD_KAPPA), -1.));

	// nsga2
	NSGA2 *n2	= new NSGA2(120, 2500);
	n2->run(dF);
	cout    << "hogehoge2" << endl << endl;
	pops	= n2->getFinalPops();
    if (n2 != NULL) {delete n2; n2 = NULL;}

	// nsga3
	// NSGA3 *n3	= new NSGA3();
	// n3->run();
	// pops	= n3->getFinalPops();
	// delete n3;

	return EXIT_SUCCESS;
}

/**
 * @fn ファイルに全個体情報を出力する
 * --- 出力形式 ---
 * # E[y1^2] E[y2y1] E[y2^2] ... E[y1^3y2^5] 目的関数1の値 目的関数2の値 ...  目的関数15の値
 * 0.0001 1.4356
 * ...
 * --- 出力形式 ---
 * @param string name ファイル名
 * @param vector<GAIndividual> &pops 個体群
 */
void Analysis::outputAllPopsIntoFile(const std::string name, const std::vector<GAIndividual> &pops)
{
    std::size_t i, ii;

	std::ofstream ofs(name);
	for (i = 0; i < pops.size(); ++i) {
		for (ii = 0; ii < pops[i].mValue.size(); ++ii) {
			ofs << pops[i].mValue[ii] << " " << std::flush;
		}
		for (ii = 0; ii < pops[i].oValue.size(); ++ii) {
			ofs << pops[i].oValue[ii] << " " << std::flush;
		}
		for (ii = 0; ii < pops[i].pValue.size(); ++ii) {
			if (ii == pops[i].pValue.size()-1) {
				ofs << pops[i].pValue[ii] << std::flush;				
			} else {
				ofs << pops[i].pValue[ii] << " " << std::flush;				
			}
		}
		ofs << std::endl;
	}
}

/**
 * main
 */
int main(int argc, char *argv[])
{
	std::string filename;

	char *ends;
	double lambda	= strtod(argv[1],&ends);
	double beta2	= strtod(argv[2],&ends);
	double alpha	= strtod(argv[3],&ends);

	std::size_t i;
	std::size_t roop = 10;	// 実験回数
	for (i = 0; i < roop; ++i) {
		std::cout << "--------------------\n" << std::endl;
		std::cout << "analysis.cpp started.\n" << std::endl;
		
		Analysis *ana	= new Analysis(lambda, beta2, alpha);

		/* GAで解く */
		std::vector<GAIndividual> pops;
		ana->GeneticAlgorithm(pops);
		filename	= "ana_gsay1pdf_" + std::to_string(i) + ".dat";
		ana->outputAllPopsIntoFile(filename, pops);

		delete ana;
		
		std::cout << "analysis.cpp has done.\n" << std::endl;
		std::cout << "--------------------\n" << std::endl;
	}
	return 0;
}
