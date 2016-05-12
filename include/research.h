#ifndef __RESEARCH_H_INCLUDE__
#define __RESEARCH_H_INCLUDE__

#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>

using namespace std;


/********** 系の係数・入力条件（不変）**********/
#define S0 (1./(2.*M_PI))	// whitenoiseのパワースペクトル
#define EPSILON 0.3			// 非線形性の強さ
#define ZETA 0.05			// 減衰定数

/********** 計算条件 **********/
#define SAMPLE_LENGTH 131072	// 131072,65536
#define NUM_OF_SAMPLES 100		// 入力の標本数
#define dt 0.01		// 時間刻み幅
#define dx 0.1		// pdfの横軸の刻み幅


class Simulation
{
private:
	double _lambda, _beta2, _alpha, _sigma;
	// 入力標本
	double _force[SAMPLE_LENGTH];
	// pdf作成時のバッファー
	double y1_buffer[SAMPLE_LENGTH][NUM_OF_SAMPLES];
	double y2_buffer[SAMPLE_LENGTH][NUM_OF_SAMPLES];
	// pdfの横軸の最小値，最大値
	double y1min, y1max, y2max, y2min;
	void _createExcitation();
	// ルンゲクッタで使う
	static double _f1(double force, double y1, double y2);
	static double _f2(double force, double y1, double y2);

public:
	Simulation(double lambda, double beta2, double alpha);
	void culcRungeKutta();
	void culcDisplacementPdf();
	void culcVelocityPdf();
	void culcDispVariance();
	void culcVelVariance();
	void exactSolutionOfGaussianWhiteNoise();
};

#endif // !__RESEARCH_H_INCLUDE_