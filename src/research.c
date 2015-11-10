#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#define PI M_PI

/********** 系の係数・入力条件（不変）**********/
#define S0 (1./(2.*PI))	// whitenoiseのパワースペクトル　0.02
#define epsilon 0.3	// 非線形性の強さ
#define zeta 0.05	// 減衰定数

/********** 入力条件（変える）**********/
//#define lambda 0.05	// パルス期待到着率
//#define beta2 10.0 // パルス振幅の分散
//#define alpha 0.1	// 入力の強さを表す:

//計算条件
#define N 131072		// 131072,65536
#define dt 0.01		// 時間刻み幅
#define dx 0.1		// pdfの横軸の刻み幅
#define SET 100		// 入力の標本数

double force[N];

// ルンゲクッタで使う
double f1(double force,double y1,double y2);
double f2(double force,double y1,double y2);


int main (int argc, char *argv[])
{
	int i, j, k, n_dx, ret1, ret2;
	double lambda, beta2, alpha;
	double wSt, pSt;
	double sigma, t, y1, y2, tmax, pdf_y1, pdf_y2, var_y1, mo4_y1, var_y2, mo4_y2, d1, d2, d3, d4;
	double DY1[4], DY2[4];
	double integral_y1 = 0., integral_y2 = 0.;
	char *ends;

	double y1_buffer[N][SET], y2_buffer[N][SET];
	double y1_pdf_buffer[3000][100], y2_pdf_buffer[3000][100];

	double y1min = 0., y1max = 0., y2max = 0., y2min = 0.;

	gsl_rng *r;
	gsl_rng *rp;

	FILE *t_force;
	FILE *t_x1;
//	FILE *x_Gpdf;

	sigma = sqrt(2.*PI*S0/dt);

	lambda = strtod(argv[1], &ends);
	beta2  = strtod(argv[2], &ends);
	alpha  = strtod(argv[3], &ends);

	wSt = sqrt(alpha);
	pSt = sqrt(1-alpha);

	t_force = fopen("t_force.dat","w");
	t_x1 	= fopen("t_x1.dat","w");
//	x_Gpdf	= fopen("y2_Gpdf.dat","w");		//ガウス分布(変位応答のpdfと比較する。)

	r  = gsl_rng_alloc (gsl_rng_default);
	rp = gsl_rng_alloc (gsl_rng_default);

	printf("The syatem is beta^2 = %lf , lambda = %lf  \n\n", beta2, lambda);
	printf("Creating a file of the pdf(.dat) ! \n\n");

	for(k=0; k<SET; k++)
	{
		gsl_rng_set (r, time (NULL)+ clock());
		gsl_rng_set (rp, time (NULL)+ clock()+1);

		// 入力を生成（ホワイトノイズ＋不規則パルス）
		for (i=0; i<N; i++)
		{
//			force[i] = (gsl_ran_gaussian(r, sigma) + gsl_ran_bernoulli(r, dt*lambda)*gsl_ran_gaussian(rp, sqrt(beta2)))/dt;
			force[i] = wSt*gsl_ran_gaussian(r, sigma) + pSt*gsl_ran_bernoulli(r, dt*lambda)*gsl_ran_gaussian(rp, sqrt(beta2))/dt;
//			force[i] = gsl_ran_gaussian(r, sigma);
//			force[i] = gsl_ran_bernoulli(r, dt*lambda)*gsl_ran_gaussian(rp, sqrt(beta2))/dt;
//			force[0] = 500;

			if(k==0)
				fprintf(t_force, "%lf %lf\n", i*dt, force[i]); //励振の記録
	}

		/*********************** Runge-Kutta **********************************/
		// 初期値
		y1 = 0.;
		y2 = 0.;
		
		for(i=0; i<N; i++)
		{
			DY1[0] = dt*f1(force[i], y1, y2);
			DY2[0] = dt*f2(force[i], y1, y2);

			DY1[1] = dt*f1(force[i], y1 + DY1[0]/2.0, y2 + DY2[0]/2.0);
			DY2[1] = dt*f2(force[i], y1 + DY1[0]/2.0, y2 + DY2[0]/2.0);

			DY1[2] = dt*f1(force[i], y1 + DY1[1]/2.0, y2 + DY2[1]/2.0);
			DY2[2] = dt*f2(force[i], y1 + DY1[1]/2.0, y2 + DY2[1]/2.0);

			DY1[3] = dt*f1(force[i], y1 + DY1[2], y2 + DY2[2]);
			DY2[3] = dt*f2(force[i], y1 + DY1[2], y2 + DY2[2]);

			y1 = y1 + (DY1[0] + 2.*DY1[1] + 2.*DY1[2] + DY1[3])/6.0;
			y2 = y2 + (DY2[0] + 2.*DY2[1] + 2.*DY2[2] + DY2[3])/6.0;

			// 変位と速度の最小値・最大値を決定
			if (y1>y1max) y1max = y1;
			if (y1<y1min) y1min = y1;
			if (y2>y2max) y2max = y2;
			if (y2<y2min) y2min = y2;

			if (i==0)
			{
				// 各入力における変位の応答を格納
				y1_buffer[i][k] = 0.;
				// 各入力における速度の応答を格納
				y2_buffer[i][k] = 0.;
			}
			else
			{
				// 各入力における変位の応答を格納
				y1_buffer[i][k] = y1;
				// 各入力における速度の応答を格納
				y2_buffer[i][k] = y2;
			}

			if(k==0)
			{
				// 応答変位の記録
				fprintf(t_x1, "%lf %lf\n", i*dt, y1);
			}		
		}
	}

	/************************* 変位のpdfを計算する ***********************************/
	FILE *y1_pdf;
	y1_pdf 	= fopen("y1_pdf.dat","w");

	for(k=0; k<SET; k++)  
	{
		for(j=0; (y1min+j*dx) <= y1max; j++)
		{
			// 幅dxに含まれる回数
			n_dx = 0;
			for(i=0; i<N; i++)
			{
				if((y1_buffer[i][k] < (y1min+(j+1)*dx)) && (y1_buffer[i][k] >= (y1min+j*dx))) n_dx++;
			}
//			printf("n_dxの中身:%lf\n",n_dx);		
			y1_pdf_buffer[j][k] = (double)n_dx/N/dx;
		}
	}

//	printf("Inside of the y1min:%lf, y1max:%lf\n",y1min, y1max);		

	for(i=0; (y1min+i*dx) <= y1max; i++)
	{
		pdf_y1 = 0.;
		for(k=0; k<SET; k++)
		{
			pdf_y1 += y1_pdf_buffer[i][k];
//			printf("y1_pdf_bufferの中身:%lf\n",y1_pdf_buffer);
		}
		pdf_y1 = (double)pdf_y1/SET;
		integral_y1 += pdf_y1*dx;
		// 応答の確率密度関数の記録
		fprintf(y1_pdf, "%lf %lf\n", (y1min+i*dx), pdf_y1);
	}

	// pdfの全積分値
	printf("integral of pdf_y1 = %lf\n",integral_y1);

	fclose(y1_pdf);

      /************************* 速度のpdfを計算する ***********************************/
//	FILE *y2_pdf;
//	y2_pdf	= fopen("y2_pdf.dat","w");
//	
//	for(k=0; k<SET; k++)  
//	{
//		for(j=0; (y2min+j*dx) <= y2max; j++ )
//		{
//			// 幅dxに含まれる回数
//			n_dx = 0;
//			for(i=0; i<N ; i++ )
//			{
//				if((y2_buffer[i][k] < (y2min+(j+1)*dx)) && (y2_buffer[i][k] >= (y2min+j*dx))) n_dx++;
//			}
//			y2_pdf_buffer[j][k] = (double)n_dx/N/dx;
//		}
//	}
//
//	for(i=0; (y2min+i*dx) <= y2max; i++)
//	{
//		pdf_y2 = 0.;
//		for(k=0; k<SET; k++)
//		{
//			pdf_y2 += y2_pdf_buffer[i][k];
//		}
//		pdf_y2 = (double)pdf_y2/SET;
//		integral_y2 += pdf_y2*dx;
//		// 応答の確率密度関数の記録
//		fprintf(y2_pdf, "%lf %lf\n", (y2min+i*dx), pdf_y2);
//	}
//
//	// pdfの全積分値
//	printf("integral of pdf_y2 = %lf\n",integral_y2);
//
//	fclose(y2_pdf);
	/************************** シミュレーション結果の変位分散を計算 *******************************/
//	FILE *y1_var;
//	y1_var  = fopen("sim_y1_var.dat", "w");
//	y1_pdf = fopen("y1_pdf.dat", "r");
//	
//	var_y1 = 0.;
//	mo4_y1 = 0.;
//	while((ret1 = fscanf(y1_pdf, "%lf %lf", &d1, &d2)) != EOF)
//	{
//		var_y1 += pow(d1,2)*d2;
//		mo4_y1 += pow(d1,4)*d2;
//	}
//	
//	fprintf(y1_var, "%lf", mo4_y1/pow(var_y1,2));
//
//	fclose(y1_pdf);
//	fclose(y1_var);
	/************************** シミュレーション結果の速度分散を計算 *******************************/
//	FILE *y2_var;
//	y2_var  = fopen("sim_y2_var.dat", "w");
//	y1_pdf = fopen("y2_pdf.dat", "r");
//
//	var_y2 = 0.;
//	while((ret2 = fscanf(y2_pdf, "%lf %lf", &d3, &d4)) != EOF)
//	{
//		var_y2 += pow(d3,2)*d4;
//	}
//	
//	fprintf(y2_var, "%lf", var_y2);
//	fclose(y2_var);
//	fclose(y2_pdf);

	/////////////////////////////////////////////////////////////
	fclose(t_force);
	fclose(t_x1);
//	fclose(x_Gpdf);

	printf ("\nIt is done !\n");
	return 0;
}

double f1(double force,double y1,double y2)
{
	return y2;
}

double f2(double force,double y1,double y2)
{
	return force - 2*zeta*y2 - y1 - epsilon*y1*y1*y1;
}

