#include "../include/common.h"

Common::Common(){}
Common::~Common(){}

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
 * @fn ２次元vectorのリサイズ
 * @param std::vector< std::vector< double > > &v ２次元vector
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