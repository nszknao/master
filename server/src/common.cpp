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
 * @param vector &x X軸情報
 * @param vector &y Y軸情報
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
	unsigned int i;
	for (i = 0; i < x.size(); ++i) {
		ofs << x[i] << " " << y[i] << std::endl;
	}
}

/**
 * @fn 2次元vectorのリサイズ
 * @param std::vector< std::vector< double > > &v 2次元vector
 * @param unsigned int s1 行数
 * @param unsigned int s2 列数
 */
void Common::resize2DemensionalVector(std::vector< std::vector< double > > &v, unsigned int s1, unsigned int s2)
{
	unsigned int i;
	v.resize(s1);
	for (i = 0; i < s1; ++i) {
		v[i].resize(s2);
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
	unsigned int i;

	for (i = 0; i < v.size(); ++i) {
		if (v[i] >= value || v[i] < -1*value)
			flg	= true;
	}

	return flg;
}
