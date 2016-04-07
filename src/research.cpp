#include "research.h"

/* 入力を生成する */
void Simulation::_createExcitation()
{
	std::cout << "Start creating excitation.\n" << std::endl;

	// カウント変数
	size_t tmp_num, tmp_lng;

	FILE *t_force;
	t_force = fopen("t_force.dat", "w");

	// ランダム変数
	gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng* rp = gsl_rng_alloc(gsl_rng_default);

	// 入力強度
	double wSt = sqrt(this->alpha);
	double pSt = sqrt(1 - this->alpha);

	for (tmp_num = 0; tmp_num < NUM_OF_SAMPLES; tmp_num++)
	{
		gsl_rng_set(r, time(NULL) + clock());
		gsl_rng_set(rp, time(NULL) + clock() + 1);

		/********** 入力を生成（ホワイトノイズ＋不規則パルス） **********/
		for (tmp_lng = 0; tmp_lng < SAMPLE_LENGTH; tmp_lng++)
		{
			this->force[tmp_lng] = wSt*gsl_ran_gaussian(r, this->sigma) + pSt*gsl_ran_bernoulli(r, dt*this->lambda)*gsl_ran_gaussian(rp, sqrt(this->beta2)) / dt;

			// 入力の記録
			if (tmp_num == 0) fprintf(t_force, "%lf %lf\n", tmp_lng*dt, this->force[tmp_lng]);
		}
	}

	fclose(t_force);
}

/* 運動方程式を4次のルンゲクッタで解く */
void Simulation::culcRungeKutta()
{
	std::cout << "Start solving equation of motion using 4th-Runge-Kutta.\n" << std::endl;

	// カウント変数
	size_t tmp_num, tmp_lng;

	FILE *t_x1;
	t_x1 = fopen("t_x1.dat", "w");

	/* 入力を生成 */
	this->_createExcitation();
	
	// ルンゲクッタの諸変数
	double DY1[4], DY2[4];

	for (tmp_num = 0; tmp_num < NUM_OF_SAMPLES; tmp_num++)
	{
		// 初期値
		double y1 = 0., y2 = 0.;

		for (tmp_lng = 0; tmp_lng < SAMPLE_LENGTH; tmp_lng++)
		{
			DY1[0] = dt*this->_f1(force[tmp_lng], y1, y2);
			DY2[0] = dt*this->_f2(force[tmp_lng], y1, y2);

			DY1[1] = dt*this->_f1(force[tmp_lng], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);
			DY2[1] = dt*this->_f2(force[tmp_lng], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);

			DY1[2] = dt*this->_f1(force[tmp_lng], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);
			DY2[2] = dt*this->_f2(force[tmp_lng], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);

			DY1[3] = dt*this->_f1(force[tmp_lng], y1 + DY1[2], y2 + DY2[2]);
			DY2[3] = dt*this->_f2(force[tmp_lng], y1 + DY1[2], y2 + DY2[2]);

			y1 = y1 + (DY1[0] + 2.*DY1[1] + 2.*DY1[2] + DY1[3]) / 6.0;
			y2 = y2 + (DY2[0] + 2.*DY2[1] + 2.*DY2[2] + DY2[3]) / 6.0;

			// 変位と速度の最小値・最大値を決定
			if (y1>y1max) y1max = y1;
			if (y1<y1min) y1min = y1;
			if (y2>y2max) y2max = y2;
			if (y2<y2min) y2min = y2;

			this->y1_buffer[tmp_lng][tmp_num] = y1;
			this->y2_buffer[tmp_lng][tmp_num] = y2;

			// 応答変位の記録
			if (tmp_num == 0) fprintf(t_x1, "%lf %lf\n", tmp_lng*dt, y1);
		}
	}

	fclose(t_x1);
}

/* 変位の確率密度関数を生成 */
void Simulation::culcDisplacementPdf()
{
	cout << "Creating a file of the displacement pdf(.dat).\n" << endl;

	// カウント変数
	size_t tmp_num, tmp_ndx, tmp_lng;

	int n_dx1;
	double y1_pdf_buffer[3000][NUM_OF_SAMPLES], pdf_y1, integral_y1 = 0.;

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
				if ((this->y1_buffer[tmp_lng][tmp_num] < (y1min + (tmp_ndx + 1)*dx)) && (this->y1_buffer[tmp_lng][tmp_num] >= (y1min + tmp_ndx*dx))) n_dx1++;
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
		}
		pdf_y1 = (double)pdf_y1 / NUM_OF_SAMPLES;
		integral_y1 += pdf_y1*dx;
		// 応答の確率密度関数の記録
		fprintf(y1_pdf, "%lf %lf\n", (y1min + tmp_ndx*dx), pdf_y1);
	}

	std::cout << "integral of pdf_y1 = " << integral_y1 << ".\n" << std::endl;

	fclose(y1_pdf);
}

/* 速度の確率密度関数を生成 */
void Simulation::culcVelocityPdf()
{
	cout << "Creating a file of the velocity pdf(.dat).\n" << endl;

	// カウント変数
	size_t tmp_num, tmp_ndx, tmp_lng;

	int n_dx;
	double y2_pdf_buffer[3000][NUM_OF_SAMPLES], pdf_y2, integral_y2 = 0.;

	FILE *y2_pdf;
	y2_pdf = fopen("y2_pdf.dat", "w");

	for (tmp_num = 0; tmp_num<NUM_OF_SAMPLES; tmp_num++)
	{
		for (tmp_ndx = 0; (y2min + tmp_ndx*dx) <= y2max; tmp_ndx++)
		{
			// 幅dxに含まれる回数
			n_dx = 0;
			for (tmp_lng = 0; tmp_lng<SAMPLE_LENGTH; tmp_lng++)
			{
				if ((this->y2_buffer[tmp_lng][tmp_num] < (y2min + (tmp_ndx + 1)*dx)) && (this->y2_buffer[tmp_lng][tmp_num] >= (y2min + tmp_ndx*dx))) n_dx++;
			}
			y2_pdf_buffer[tmp_ndx][tmp_num] = (double)n_dx / SAMPLE_LENGTH / dx;
		}
	}

	for (tmp_ndx = 0; (y2min + tmp_ndx*dx) <= y2max; tmp_ndx++)
	{
		pdf_y2 = 0.;
		for (tmp_num = 0; tmp_num<NUM_OF_SAMPLES; tmp_num++)
		{
			pdf_y2 += y2_pdf_buffer[tmp_ndx][tmp_num];
		}
		pdf_y2 = (double)pdf_y2 / NUM_OF_SAMPLES;
		integral_y2 += pdf_y2*dx;
		// 応答の確率密度関数の記録
		fprintf(y2_pdf, "%lf %lf\n", (y2min + tmp_ndx*dx), pdf_y2);
	}

	// pdfの全積分値
	std::cout << "integral of pdf_y2 = " << integral_y2 << ".\n" << std::endl;

	fclose(y2_pdf);
}

/* 変位分散を計算 */
void Simulation::culcDispVariance()
{
	cout << "Creating a file of the displacement varience(.dat).\n" << endl;

	int scan_file;
	double var_y1 = 0., mo4_y1 = 0., first_row, second_row;
	
	FILE *y1_var, *y1_pdf;
	y1_var = fopen("sim_y1_var.dat", "w");
	y1_pdf = fopen("y1_pdf.dat", "r");

	while ((scan_file = fscanf(y1_pdf, "%lf %lf", &first_row, &second_row)) != EOF)
	{
		var_y1 += pow(first_row, 2)*second_row;
		mo4_y1 += pow(first_row, 4)*second_row;
	}

	fprintf(y1_var, "%lf", mo4_y1 / pow(var_y1, 2));

	fclose(y1_pdf);
	fclose(y1_var);
}

/* 速度分散を計算 */
void Simulation::culcVelVariance()
{
	cout << "Creating a file of the velocity varience(.dat).\n" << endl;

	int scan_file;
	double var_y2 = 0., first_row, second_row;

	FILE *y2_var, *y2_pdf;
	y2_var = fopen("sim_y2_var.dat", "w");
	y2_pdf = fopen("y2_pdf.dat", "r");

	while ((scan_file = fscanf(y2_pdf, "%lf %lf", &first_row, &second_row)) != EOF)
	{
		var_y2 += pow(first_row, 2)*second_row;
	}

	fprintf(y2_var, "%lf", var_y2);

	fclose(y2_var);
	fclose(y2_pdf);
}

/* ガウス性ホワイトノイズを受ける系の厳密解 */
void Simulation::exactSolutionOfGaussianWhiteNoise()
{
	cout << "Culculate exact solution excited Gaussian white noise." << endl;

	// カウント変数
	size_t tmp;

	FILE *x_Gpdf;
	x_Gpdf = fopen("y1_Gpdf.dat", "w");

	// 第二修正ベッセル関数
	double K_nu;
	double bessel_p, bessel_q;
	double bessel_C, bessel_arg;

	// 確率密度関数の厳密解epsilon
	double exact_gauss_pdf;
	double integral_gauss_pdf;

	bessel_p = 2.*ZETA;
	bessel_q = ZETA*EPSILON;
	bessel_arg = pow(bessel_p, 2) / (8.*bessel_q);

	K_nu = gsl_sf_bessel_Knu(0.25, bessel_arg);
	bessel_C = 2.*sqrt(bessel_q / bessel_p)*exp(-bessel_arg) / K_nu;

	for (tmp = 0; (y1min + tmp*dx) <= y1max; tmp++)
	{
		exact_gauss_pdf = bessel_C*exp(-bessel_p*pow(y1min + tmp*dx, 2) - bessel_q*pow(y1min + tmp*dx, 4));
		integral_gauss_pdf += exact_gauss_pdf*dx;
		// 厳密解の記録
		fprintf(x_Gpdf, "%lf %lf\n", (y1min + tmp*dx), exact_gauss_pdf);
	}

	std::cout << "integral of exact_gauss_pdf = " << integral_gauss_pdf << ".\n" << std::endl;

	fclose(x_Gpdf);
}

double Simulation::_f1(double force, double y1, double y2)
{
	return y2;
}

double Simulation::_f2(double force, double y1, double y2)
{
	return force - 2*ZETA*y2 - y1 - EPSILON*y1*y1*y1;
}

int main (int argc, char *argv[]) {

	Simulation *sim;
	sim = new Simulation;

	// 入力パラメータ
	char *ends;
	sim->lambda = strtod(argv[1], &ends);
	sim->beta2 = strtod(argv[2], &ends);
	sim->alpha = strtod(argv[3], &ends);
	sim->sigma = sqrt(2.*M_PI*S0 / dt);

	cout << "--------------------------\n" << endl;
	cout << "research.cpp started.\n" << endl;

	sim->culcRungeKutta();
	sim->culcDisplacementPdf();
	sim->exactSolutionOfGaussianWhiteNoise();

	cout << "research.cpp has done.\n" << endl;
	cout << "--------------------------\n" << endl;
	return 0;
}
