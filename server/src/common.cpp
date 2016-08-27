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
 * @fn 特定のvectorを基準にしてソートする
 * @param vector< vector<double> > &ref 基準となる配列
 * @param vector< vector<double> > &target ソートしたい配列
 * @param int key ソートしたいキー
 */
void Common::sortBasedOnParticularArray(const std::vector< std::vector<double> > &ref, std::vector< std::vector<double> > &target, int key)
{
	if (ref.size() != target.size()) {
		std::cout << "Error: Different size of vector.";
		exit(1);
	}

	unsigned int i;
	struct data_t {
		int index;
		double value;
	};

	std::vector<data_t> refValue(ref.size());
	for (i = 0; i < refValue.size(); ++i) {
		refValue[i].index	= i;
		refValue[i].value	= ref[i][key];
	}
	// for (i = 0; i < refValue.size(); ++i) {
	// 	std::cout << "before [ " << refValue[i].index << " : " << refValue[i].value << " ]" << std::endl;
	// }
	sort(refValue.begin(), refValue.end(), [](const data_t &a, const data_t &b){
		return (a.value == b.value) ? (a.index < b.index) : (a.value < b.value);
	});
	// for (i = 0; i < refValue.size(); ++i) {
	// 	std::cout << "after [ " << refValue[i].index << " : " << refValue[i].value << " ]" << std::endl;
	// }
	// 対象の配列をソート
	std::vector< std::vector<double> > tmp(target.size());
	for (i = 0; i < target.size(); ++i) {
		tmp[i]	= target[refValue[i].index];
	}
	target	= tmp;
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