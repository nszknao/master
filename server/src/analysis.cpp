#include <analysis.h>
#include <expfit.h>
#include <nsga2.h>
#include <nsga3.h>
#include <common.h>
#include <ga_individual.h>

const double Analysis::GGD_KAPPA = 2.;  // 1.:ラプラス分布，2.:ガウス分布，∞:一様分布

Analysis::Analysis(double arg_l, double arg_b, double arg_a)
:_lambda(arg_l),
_beta2(arg_b),
_alpha(arg_a)
{}

Analysis::~Analysis() {}

/**
 * @fn NSGA2でモーメント方程式を解く
 * @param vector<double> &pValue 計算結果のパラメータ値を保存
 * @param vector<double> &oValue 計算結果の目的関数値を保存
 */
std::vector<GAIndividual> Analysis::GeneticAlgorithm(std::size_t pop, std::size_t gen)
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
    
    // 目的関数を選ぶ
    std::vector<std::size_t> selectedObj{1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14};

    // 目的関数の重み計算用NSGA2
//    MomentEq meq1;
//    meq1.setPrmdG(dF);
//    meq1.setObjList(selectedObj);
//    NSGA2 *n2_  = new NSGA2(1000, 0);
//    n2_->run(&meq1);
//    std::vector<GAIndividual> pops_ = n2_->getFinalPops();
//    if (n2_ != NULL) {delete n2_; n2_ = NULL;}
//    std::vector< std::vector<double> > objWeight = this->_getNormalizeObjectList(pops_);
   
    MomentEq meq2;
    meq2.setPrmdG(dF);
    meq2.setObjList(selectedObj);
//    meq2.setObjWeight(objWeight);

    // NSGA2 or NSGA3
//    NSGA2 *nsga = new NSGA2(pop, gen);
//    NSGA3 *nsga = new NSGA3();
//    jMetalCpp_ssNSGAII *nsga = new jMetalCpp_ssNSGAII(pop, gen);
//    jMetalCpp_NSGAII *nsga = new jMetalCpp_NSGAII(pop, gen);
    jMetalCpp_GDE3 *nsga = new jMetalCpp_GDE3(pop, gen);

    nsga->run(&meq2);
    std::vector<GAIndividual> pops = nsga->getFinalPops();
    if (nsga != NULL) {delete nsga; nsga = NULL;}

    return pops;
}

/**
 * @fn 個体群から目的関数の重み用の標準偏差リストを作成する
 * @param vector<GAIndividual> &pops 個体群
 * @return vector< vector<double> > list 目的関数の正規化用リスト（0: 平均，1: 標準偏差）
 */
std::vector< std::vector<double> >
Analysis::_getNormalizeObjectList(const std::vector<GAIndividual> &pops)
{
    std::size_t i, ii;
    double sumObj, sumSquareObj;

    std::vector< std::vector<double> > list(2, std::vector<double>(Common::NUM_OF_MOMENTEQ, 0));
    
    for (i = 0; i < Common::NUM_OF_MOMENTEQ; ++i) {
        // 平均
        sumObj = 0.;
        for (ii = 0; ii < pops.size(); ++ii) {
            sumObj += pops[ii].oValue[i];
        }
        list[0][i] = sumObj / pops.size();
        // 標準偏差
        sumSquareObj = 0.;
        for (ii = 0; ii < pops.size(); ++ii) {
            sumSquareObj += pow((pops[ii].oValue[i] - list[0][i]), 2);
        }
        list[1][i] = sqrt(sumSquareObj / (pops.size()-1.));
    }
    return list;
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
    double lambda   = strtod(argv[1],&ends);
    double beta2    = strtod(argv[2],&ends);
    double alpha    = strtod(argv[3],&ends);

    std::size_t i;
    std::size_t roop = 5;  // 実験回数
    for (i = 0; i < roop; ++i) {
        std::cout << "--------------------" << std::endl;
        std::cout << "analysis.cpp started.\n" << std::endl;
        
        try {
            Analysis *ana   = new Analysis(lambda, beta2, alpha);

            /* GAで解く */
            std::vector<GAIndividual> pops;
            if (i == 0) {
                //pops = ana->GeneticAlgorithm(100, 25000);
                pops = ana->GeneticAlgorithm(100, 230);
            }
            else if (i == 1) {
                //pops = ana->GeneticAlgorithm(200, 25000);
                pops = ana->GeneticAlgorithm(100, 240);
            }
            else if (i == 2) {
                //pops = ana->GeneticAlgorithm(300, 25000);
                pops = ana->GeneticAlgorithm(100, 250);
            }
            else if (i == 3) {
                //pops = ana->GeneticAlgorithm(400, 25000);
                pops = ana->GeneticAlgorithm(100, 260);
            }
            else if (i == 4) {
                //pops = ana->GeneticAlgorithm(500, 25000);
                pops = ana->GeneticAlgorithm(100, 270);
            }
            filename    = "ana_gsay1pdf_" + std::to_string(i) + ".dat";
            ana->outputAllPopsIntoFile(filename, pops);
    
            delete ana;
        } catch(std::exception e) {
            std::cout << e.what() << std::endl;
        } catch(SharkException e) {
            std::cout << e.what() << std::endl;
        }
                
        std::cout << "analysis.cpp has done.\n" << std::endl;
        std::cout << "--------------------\n" << std::endl;
    }
    return 0;
}
