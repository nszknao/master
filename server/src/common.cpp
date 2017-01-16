#include "../include/common.h"

Common::Common(){}
Common::~Common(){}

/********** 系の係数・入力条件**********/
const double Common::S0 = 1./(2.*M_PI);
const double Common::EPSILON = 0.3;
const double Common::ZETA = 0.05;

/********** 解析条件 **********/
const std::size_t Common::NUM_OF_MOMENTEQ = 15;
const std::size_t Common::NUM_OF_MOMENT = 21;
const std::size_t Common::NUM_OF_PARAM = 10;

/**
 * @fn ファイルに出力する
 * @param string name ファイル名
 * @param vector<double> &x X軸情報
 * @param vector<double> &y Y軸情報
 */
void Common::outputIntoFile(const std::string name, const std::vector<double> &x, const std::vector<double> &y)
{
    std::cout << "Creating a file.\n" << std::endl;

    // 値の数をチェック
    if (x.size() != y.size()) {
        std::cout << "Error: Do not match vector size." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ofstream ofs(name);
    std::size_t i;
    for (i = 0; i < x.size(); ++i) {
        ofs << x[i] << " " << y[i] << std::endl;
    }
}

/**
 * @fn ファイルに出力する
 * @param string name ファイル名
 * @param vector<double> &x X軸情報
 * @param vector<double> &y Y軸情報
 * @param vector<double> &z Z軸情報
 */
void Common::output3DIntoFile(const std::string name, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z)
{
    std::cout << "Creating a file.\n" << std::endl;

    // 値の数をチェック
    if ((x.size() != y.size()) || (pow(y.size(),2) != z.size()) || (z.size() != pow(x.size(),2))) {
        std::cout << "Error: Do not match vector size." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ofstream ofs(name);
    std::size_t i, ii;
    // 1行目
    ofs << "X" << ",";
    for (i = 0; i < x.size(); ++i) {
        if (i == x.size()-1) {
            ofs << x[i] << std::endl;
        } else {
            ofs << x[i] << ",";
        }
    }
    // 2行目以降
    for (i = 0; i < y.size(); ++i) {
        ofs << y[i] << ",";
        for (ii = 0; ii < x.size(); ++ii) {
            if (ii == x.size()-1) {
                ofs << z[i*x.size() + ii] << std::endl;
            } else {
                ofs << z[i*x.size() + ii] << ",";
            }
        }
    }
}

/**
 * @fn 2次元vectorのリサイズ
 * @param vector< vector< double > > &v 2次元vector
 * @param size_t s1 行数
 * @param size_t s2 列数
 */
void Common::resize2DemensionalVector(std::vector< std::vector< double > > &v, std::size_t s1, std::size_t s2)
{
    std::size_t i;
    v.resize(s1);
    for (i = 0; i < s1; ++i) {
        v[i].resize(s2);
    }
}

/**
 * @fn 3次元vectorのリサイズ
 * @param vector< vector< vector<double> > > &v 3次元vector
 * @param size_t s1 行数
 * @param size_t s2 列数
 * @param size_t s3 何か
 */
void Common::resize3DemensionalVector(std::vector< std::vector< std::vector<double> > > &v, std::size_t s1, std::size_t s2, std::size_t s3)
{
    std::size_t i, ii;
    v.resize(s1);
    for (i = 0; i < s1; ++i) {
        v[i].resize(s2);
        for (ii = 0; ii < s2; ++ii) {
            v[i][ii].resize(s3);
        }
    }
}

/**
 * @fn 指定した配列の要素が閾値以内に収まっているか判定
 * @param vector<double> &v 対象のvector
 * @param double value 閾値
 */
bool Common::isOverSpecifyValue(const std::vector<double> &v, double value)
{
    bool flg = false;
    std::size_t i;

    for (i = 0; i < v.size(); ++i) {
        if (v[i] >= value || v[i] < -1*value)
            flg = true;
    }

    return flg;
}
