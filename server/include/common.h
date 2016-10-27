#ifndef __COMMON_H_INCLUDE__
#define __COMMON_H_INCLUDE__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

class Common
{
private:
public:
	/********** 系の係数・入力条件**********/
	static const double S0;
	static const double EPSILON;
	static const double ZETA;

	/********* 解析条件 *********/
	static const std::size_t NUM_OF_MOMENTEQ;
	static const std::size_t NUM_OF_MOMENT;
	static const std::size_t NUM_OF_PARAM;

	Common();
	~Common();
	static void outputIntoFile(const std::string, const std::vector<double> &, const std::vector<double> &);
	static void resize2DemensionalVector(std::vector< std::vector< double > > &, unsigned int, unsigned int);
	static bool isOverSpecifyValue(const std::vector<double> &, double);
};

#endif // !__COMMON_H_INCLUDE_