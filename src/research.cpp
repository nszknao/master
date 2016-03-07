#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>


/********** 系の係数・入力条件（不変）**********/
#define S0 (1./(2.*M_PI))	// whitenoiseのパワースペクトル
#define epsilon 0.3			// 非線形性の強さ
#define zeta 0.05			// 減衰定数

/********** 計算条件 **********/
#define SAMPLE_LENGTH 131072	// 131072,65536
#define NUM_OF_SAMPLES 100		// 入力の標本数
#define dt 0.01		// 時間刻み幅
#define dx 0.1		// pdfの横軸の刻み幅

// ルンゲクッタで使う
double f1(double force, double y1, double y2);
double f2(double force, double y1, double y2);


int main (int argc, char *argv[]) {

	// 入力パラメータ
	char *ends;
	double lambda	= strtod(argv[1], &ends);
	double beta2	= strtod(argv[2], &ends);
	double alpha	= strtod(argv[3], &ends);
	double sigma	= sqrt(2.*M_PI*S0 / dt);

	// カウント変数
	size_t tmp, tmp_ndx, tmp_num, tmp_lng;
	// 入力標本
	double force[SAMPLE_LENGTH];

	gsl_rng *r;
	gsl_rng *rp;

	FILE *t_force;
	FILE *t_x1;

	t_force = fopen("t_force.dat", "w");
	t_x1 = fopen("t_x1.dat", "w");

	r = gsl_rng_alloc(gsl_rng_default);
	rp = gsl_rng_alloc(gsl_rng_default);
    
    cout << "The syatem is beta^2 = " << beta2 << ", lambda = " << lambda << "\n\n" << endl;
    cout << "Creating a file of the pdf(.dat)!\n\n" endl;
	
	// 入力強度
	double wSt = sqrt(alpha);
	double pSt = sqrt(1 - alpha);

	// ルンゲクッタで使う
	double DY1[4], DY2[4];
	double y1_buffer[SAMPLE_LENGTH][NUM_OF_SAMPLES] = {{0.}}, y2_buffer[SAMPLE_LENGTH][NUM_OF_SAMPLES] = {{0.}};
	double y1min = 0., y1max = 0., y2max = 0., y2min = 0.;

	for (tmp_num = 0; tmp_num < NUM_OF_SAMPLES; tmp_num++)
	{
		gsl_rng_set(r, time(NULL) + clock());
		gsl_rng_set(rp, time(NULL) + clock() + 1);

		/********** 入力を生成（ホワイトノイズ＋不規則パルス） **********/
		for (tmp_lng = 0; tmp_lng < SAMPLE_LENGTH; tmp_lng++)
		{
			force[tmp_lng] = wSt*gsl_ran_gaussian(r, sigma) + pSt*gsl_ran_bernoulli(r, dt*lambda)*gsl_ran_gaussian(rp, sqrt(beta2)) / dt;
			// force[tmp_lng] = gsl_ran_gaussian(r, sigma);
			// force[tmp_lng] = gsl_ran_bernoulli(r, dt*lambda)*gsl_ran_gaussian(rp, sqrt(beta2))/dt;

			// 入力の記録
			if (tmp_num == 0) fprintf(t_force, "%lf %lf\n", tmp_lng*dt, force[tmp_lng]);
		}

		/*********************** Runge-Kutta **********************************/
		// 初期値
		double y1 = 0.;
		double y2 = 0.;

		for (tmp_lng = 0; tmp_lng < SAMPLE_LENGTH; tmp_lng++)
		{
			DY1[0] = dt*f1(force[tmp_lng], y1, y2);
			DY2[0] = dt*f2(force[tmp_lng], y1, y2);

			DY1[1] = dt*f1(force[tmp_lng], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);
			DY2[1] = dt*f2(force[tmp_lng], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);

			DY1[2] = dt*f1(force[tmp_lng], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);
			DY2[2] = dt*f2(force[tmp_lng], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);

			DY1[3] = dt*f1(force[tmp_lng], y1 + DY1[2], y2 + DY2[2]);
			DY2[3] = dt*f2(force[tmp_lng], y1 + DY1[2], y2 + DY2[2]);

			y1 = y1 + (DY1[0] + 2.*DY1[1] + 2.*DY1[2] + DY1[3]) / 6.0;
			y2 = y2 + (DY2[0] + 2.*DY2[1] + 2.*DY2[2] + DY2[3]) / 6.0;

			// 変位と速度の最小値・最大値を決定
			if (y1>y1max) y1max = y1;
			if (y1<y1min) y1min = y1;
			if (y2>y2max) y2max = y2;
			if (y2<y2min) y2min = y2;

			y1_buffer[tmp_lng][tmp_num] = y1;
			y2_buffer[tmp_lng][tmp_num] = y2;

			// 応答変位の記録
			if (tmp_num == 0) fprintf(t_x1, "%lf %lf\n", tmp_lng*dt, y1);
		}
	}

	/************************* 変位のpdfを計算する ***********************************/
	int n_dx1;
	double y1_pdf_buffer[3000][100], pdf_y1, integral_y1 = 0.;
	FILE *y1_pdf;
	y1_pdf = fopen("y1_pdf.dat", "w");

	for (tmp_num = 0; tmp_num < NUM_OF_SAMPLES; tmp_num++)
	{
		for (tmp_ndx = 0; (y1min + tmp_ndx*dx) <= y1max; tmp_ndx++)
		{
			// 幅dxに含まれる回数
			n_dx1 = 0;
			for (tmp_lng = 0; tmp_lng < SAMPLE_LENGTH; tmp_lng++)
			{
				if ((y1_buffer[tmp_lng][tmp_num] < (y1min + (tmp_ndx + 1)*dx)) && (y1_buffer[tmp_lng][tmp_num] >= (y1min + tmp_ndx*dx))) n_dx1++;
			}
			y1_pdf_buffer[tmp_ndx][tmp_num] = (double)n_dx1 / SAMPLE_LENGTH / dx;
		}
	}

	for (tmp_ndx = 0; (y1min + tmp_ndx*dx) <= y1max; tmp_ndx++)
	{
		pdf_y1 = 0.;
		for (tmp_num = 0; tmp_num < NUM_OF_SAMPLES; tmp_num++)
		{
			pdf_y1 += y1_pdf_buffer[tmp_ndx][tmp_num];
			// printf("y1_pdf_bufferの中身:%lf\n",y1_pdf_buffer);
		}
		pdf_y1 = (double)pdf_y1 / NUM_OF_SAMPLES;
		integral_y1 += pdf_y1*dx;
		// 応答の確率密度関数の記録
		fprintf(y1_pdf, "%lf %lf\n", (y1min + tmp_ndx*dx), pdf_y1);
	}

	printf("integral of pdf_y1 = %lf\n", integral_y1);

	fclose(y1_pdf);

	/************************* 速度のpdfを計算する ***********************************/
	//	int n_dx2, pdf_y2, integral_y2 = 0;
	//	double y2_pdf_buffer[3000][100];
	//	FILE *y2_pdf;
	//	y2_pdf	= fopen("y2_pdf.dat","w");
	//	
	//	for(tmp_num=0; tmp_num<NUM_OF_SAMPLES; tmp_num++)  
	//	{
	//		for(tmp_ndx=0; (y2min+tmp_ndx*dx) <= y2max; tmp_ndx++ )
	//		{
	//			// 幅dxに含まれる回数
	//			n_dx2 = 0;
	//			for(tmp_lng=0; tmp_lng<SAMPLE_LENGTH ; tmp_lng++ )
	//			{
	//				if((y2_buffer[tmp_lng][tmp_num] < (y2min+(tmp_ndx+1)*dx)) && (y2_buffer[tmp_lng][tmp_num] >= (y2min+tmp_ndx*dx))) n_dx++;
	//			}
	//			y2_pdf_buffer[tmp_ndx][tmp_num] = (double)n_dx2/SAMPLE_LENGTH/dx;
	//		}
	//	}
	//
	//	for(tmp_ndx=0; (y2min+tmp_ndx*dx) <= y2max; tmp_ndx++)
	//	{
	//		pdf_y2 = 0.;
	//		for(tmp_num=0; tmp_num<NUM_OF_SAMPLES; tmp_num++)
	//		{
	//			pdf_y2 += y2_pdf_buffer[tmp_ndx][tmp_num];
	//		}
	//		pdf_y2 = (double)pdf_y2/NUM_OF_SAMPLES;
	//		integral_y2 += pdf_y2*dx;
	//		// 応答の確率密度関数の記録
	//		fprintf(y2_pdf, "%lf %lf\n", (y2min+tmp_ndx*dx), pdf_y2);
	//	}
	//
	//	// pdfの全積分値
	//	printf("integral of pdf_y2 = %lf\n",integral_y2);
	//
	//	fclose(y2_pdf);

	/************************** シミュレーション結果の変位分散を計算 *******************************/
	//	int scan_file, var_y1, mo4_y1, first_row, second_row;
	//	FILE *y1_var;
	//	y1_var  = fopen("sim_y1_var.dat", "w");
	//	y1_pdf = fopen("y1_pdf.dat", "r");
	//	
	//	var_y1 = 0.;
	//	mo4_y1 = 0.;
	//	while((scan_file = fscanf(y1_pdf, "%lf %lf", &first_row, &second_row)) != EOF)
	//	{
	//		var_y1 += pow(first_row,2)*second_row;
	//		mo4_y1 += pow(first_row,4)*second_row;
	//	}
	//	
	//	fprintf(y1_var, "%lf", mo4_y1/pow(var_y1,2));
	//
	//	fclose(y1_pdf);
	//	fclose(y1_var);

	/************************** シミュレーション結果の速度分散を計算 *******************************/
	//	int scan_file, var_y2, first_row, second_row;
	//	int var_y2;
	//	FILE *y2_var;
	//	y2_var  = fopen("sim_y2_var.dat", "w");
	//	y1_pdf = fopen("y2_pdf.dat", "r");
	//
	//	var_y2 = 0.;
	//	while((scan_file = fscanf(y2_pdf, "%lf %lf", &first_row, &second_row)) != EOF)
	//	{
	//		var_y2 += pow(first_row,2)*second_row;
	//	}
	//	
	//	fprintf(y2_var, "%lf", var_y2);
	//	fclose(y2_var);
	//	fclose(y2_pdf);

	/*************** ガウス性ホワイトノイズの厳密解 ***************/
	//	FILE *x_Gpdf;
	//	x_Gpdf	= fopen("y1_Gpdf.dat","w");
	//	// 第二修正ベッセル関数
	//	double K_nu;
	//	double bessel_p, bessel_q;
	//	double bessel_C, bessel_arg;
	//	// 確率密度関数の厳密解
	//	double exact_gauss_pdf;
	//	double integral_gauss_pdf;
	//
	//	bessel_p	= 2*zeta;
	//	bessel_q	= zeta*epsilon;
	//	bessel_arg	= pow(bessel_p,2)/(8*bessel_q);
	//
	//	K_nu	= gsl_sf_bessel_Knu(0.25, bessel_arg);
	//	bessel_C	= sqrt(bessel_q/bessel_p)*exp(-bessel_arg)/K_nu;
	//	
	//	for(i=0; (y1min+i*dx) <= y1max; i++)
	//	{
	//		exact_gauss_pdf = 0.;
	//		exact_gauss_pdf	= 2*bessel_C*exp(-bessel_p*pow(y1min+i*dx,2) - bessel_q*pow(y1min+i*dx,4));
	//		integral_gauss_pdf += exact_gauss_pdf*dx;
	//		// 厳密解の記録
	//		fprintf(x_Gpdf, "%lf %lf\n", (y1min+i*dx), exact_gauss_pdf);
	//	}
	//	printf("integral of exact_gauss_pdf = %f\n", integral_gauss_pdf);
	//
	//	fclose(x_Gpdf);

	fclose(t_force);
	fclose(t_x1);

	printf("\nresearch.cpp has done !\n");
	return 0;
}

double f1(double force, double y1, double y2)
{
	return y2;
}

double f2(double force, double y1, double y2)
{
	return force - 2 * zeta*y2 - y1 - epsilon*y1*y1*y1;
}

