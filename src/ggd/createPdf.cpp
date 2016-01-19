/*
 * @brief 一般化ガウス分布（generalized Gauss distribution）
 * を生成してファイルに書き出す．
 * 
 * ****パラメータ****
 * main関数の引数・・・		beta2
 * 一般化ガウス分布の係数・・・	GGD_KAPPA
 * 出力するファイル名・・・		file_name
 */
#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>


/********** 系の係数・入力条件（不変）**********/
#define GGD_KAPPA 1.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布

/********** 計算条件 **********/
#define SAMPLE_LENGTH 131072	// 131072,65536
#define dt 0.01		// 時間刻み幅

int main (int argc, char *argv[]) {

	// 入力パラメータ
	char *ends;
	double beta2	= strtod(argv[1], &ends);

	// カウント変数
	size_t tmp_dx;

	// 入力標本
	double force[SAMPLE_LENGTH];

	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a		= sqrt(gsl_sf_gamma(1./GGD_KAPPA)*pow(gsl_sf_gamma(3./GGD_KAPPA), -1.)*beta2);

	double x;
	/********** 入力を生成（ホワイトノイズ＋不規則パルス） **********/
	// 出力ファイル名
	string file_name	= "laplace_t_force.dat";

	FILE *t_force;
	t_force = fopen(file_name, "w");
	#define x_min -30

	for (tmp_dx = 0; tmp_dx*dt < 60; tmp_dx++)
	{
		x	= (double)x_min + tmp_dx*dt;
		force[tmp_dx]	= gsl_ran_exppow_pdf(x, ggd_a, GGD_KAPPA);
		// 入力の記録
		fprintf(t_force, "%lf %lf\n", x, force[tmp_dx]);
	}

	fclose(t_force);

	return 0;
}