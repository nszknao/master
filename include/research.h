#ifndef __RESEARCH_H_INCLUDE__
#define __RESEARCH_H_INCLUDE__

#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>

#include "common.h"

using namespace std;


/********** 系の係数・入力条件（不変）**********/
#define S0 (1./(2.*M_PI))	// whitenoiseのパワースペクトル
#define EPSILON 0.3			// 非線形性の強さ
#define ZETA 0.05			// 減衰定数

/********** 計算条件 **********/
#define SAMPLE_LENGTH 131072	// 131072,65536
#define NUM_OF_SAMPLES 200		// 入力の標本数
#define dt 0.01		// 時間刻み幅
#define dx 0.1		// pdfの横軸の刻み幅


class Simulation
{
private:
	double _lambda, _beta2, _alpha, _sigma;
	// pdfの横軸の最小値，最大値
	double _y1min, _y1max, _y2max, _y2min;
	// ルンゲクッタで使う
	static double _f1(double force, double y1, double y2);
	static double _f2(double force, double y1, double y2);

public:
	Simulation(double lambda, double beta2, double alpha);
	void createExcitation(std::vector<double> &, std::vector<double> &);
	void culcRungeKutta(std::vector<double> &, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &, const std::vector<double> &);
	void createDispPdf(const std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double> &);
	void createVelPdf(const std::vector< std::vector<double> > &, std::vector<double> &, std::vector<double> &);
	void exactSolutionOfGaussianWhiteNoise();
};

#endif // !__RESEARCH_H_INCLUDE_