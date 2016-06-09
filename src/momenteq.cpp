/*
	モーメント方程式

	パラメータベクトルx
	std::vector<double> x = {a, μ1, μ2, σ11, σ12, σ21, σ22, k1, k2, k3}
*/

#include "../include/momentEq.h"

using namespace std;

void momentEqFunc::momentEqFunc(const std::vector<double>& x, paramdata *params)
{
	int NUM_OF_MOMENT_EQUATION	= paramData->n;
	double *y				= paramData->y;
	double zeta				= paramData->zeta;
	double epi				= paramData->epsilon;
	double *dG				= paramData->dG;

	double cf_moment_eq[] =						// モーメント方程式 係数行列 15x21
	{
		0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		-1,-2*zeta,1,-1*epi,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,-2,-4*zeta,0,-2*epi,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	       	
		0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,-1,-2*zeta,3,0,0,-1*epi,0,0,0,0,0,0,0,0,0,0,0,0,
		dG[1],0,0,0,-2,-4*zeta,2,0,0,-2*epi,0,0,0,0,0,0,0,0,0,0,0, 
		0,3*dG[1],0,0,0,-3,-6*zeta,1,0,0,-3*epi,0,0,0,0,0,0,0,0,0,0,
		0,0,6*dG[1],0,0,0,-4,-8*zeta,0,0,0,-4*epi,0,0,0,0,0,0,0,0,0,

		0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,-1,-2*zeta,5,0,0,0,0,-1*epi,0,0,0,0,0,
		0,0,0,dG[1],0,0,0,0,0,-2,-4*zeta,4,0,0,0,0,-2*epi,0,0,0,0,
		0,0,0,0,3*dG[1],0,0,0,0,0,-3,-6*zeta,3,0,0,0,0,-3*epi,0,0,0,
		dG[3],0,0,0,0,6*dG[1],0,0,0,0,0,-4,-8*zeta,2,0,0,0,0,-4*epi,0,0,
		0,5*dG[3],0,0,0,0,10*dG[1],0,0,0,0,0,-5,-10*zeta,1,0,0,0,0,-5*epi,0,
		0,0,15*dG[3],0,0,0,0,15*dG[1],0,0,0,0,0,-6,-12*zeta,0,0,0,0,0,-6*epi
	};
	
	// モーメント方程式を解く際の補正係数
	double keisu[15];
	keisu[0] = 5.;
	keisu[1] = 24. / 5.;
	keisu[2] = 24. / 6.;
	keisu[3] = 24. / 16.;
	keisu[4] = 24. / 19.;
	keisu[5] = 1.;
	keisu[6] = 24. / 32.;
	keisu[7] = 24. / 36.;
	keisu[8] = 24. / 64.;
	keisu[9] = 24. / 79.;
	keisu[10] = 24. / 98.;
	keisu[11] = 24. / 117.;
	keisu[12] = 24. / 144;
	keisu[13] = 24. / 183.;
	keisu[14] = 24. / 216;

	// カウント変数
	size_t tmp;
	
	// ゆとり作戦
	double pa0 = x[0], pa1 = x[1], pa2 = x[2], pa3 = x[3], pa4 = x[4], pa5 = x[5], pa6 = x[6], pa7 = x[7], pa8 = x[8], pa9 = x[9];

	/******************** モーメント（開始） ********************/
	double Eg[21];

	// ２次モーメント
	Eg[0] = (1 - pa0)*(pow(pa3,2) + pow(pa1,2)) + pa0*pow(pa5,2);														// y_1^2
	Eg[1] = (1 - pa0)/2*((pa7 + pa1*pa2) + (pa9 + pa1*pa2)) + pa0*pa8;													// y_1*y_2
	Eg[2] = (1 - pa0)*(pow(pa4,2) + pow(pa2,2)) + pa0*pow(pa6,2);														// y_2^2
	
	// ４次モーメント
	Eg[3] = (1 - pa0)*(3*pow(pa3,4) + 6*pow(pa1,2)*pow(pa3,2) + pow(pa1,4)) + 3*pa0*pow(pa5,4);										// y_1^4
	Eg[4] = (1 - pa0)/2*((3*pow(pa3,2)*(pa7 + pa1*pa2) + pow(pa1,2)*(3*pa7 + pa1*pa2))
		+ (3*pow(pa3,2)*(pa9 + pa1*pa2) + pow(pa1,2)*(3*pa9 + pa1*pa2)))
		+ 3*pa0*pow(pa5,2)*pa8;																		// y_1^3*y_2
	Eg[5] = (1 - pa0)/2*((pow(pa3,2)*pow(pa4,2) + pow(pa1,2)*pow(pa4,2) + pow(pa2,2)*pow(pa3,2) + 2*pow(pa7,2) + 4*pa7*pa1*pa2 + pow(pa1,2)*pow(pa2,2))
		+ (pow(pa3,2)*pow(pa4,2) + pow(pa1,2)*pow(pa4,2) + pow(pa2,2)*pow(pa3,2) + 2*pow(pa9,2) + 4*pa9*pa1*pa2 + pow(pa1,2)*pow(pa2,2)))
		+ pa0*(pow(pa5,2)*pow(pa6,2) + 2*pow(pa8,2));															// y_1^2*y_2^2
	Eg[6] = (1 - pa0)/2*((3*pow(pa4,2)*(pa7 + pa1*pa2) + pow(pa2,2)*(3*pa7 + pa1*pa2))
		+ (3*pow(pa4,2)*(pa9 + pa1*pa2) + pow(pa2,2)*(3*pa9 + pa1*pa2)))
		+ 3*pa0*pow(pa6,2)*pa8;																		// y_1*y_2^3
	Eg[7] = (1 - pa0)*(3*pow(pa4,4) + 6*pow(pa2,2)*pow(pa4,2) + pow(pa2,4)) + 3*pa0*pow(pa6,4);										// y_2^4
	
	// ６次モーメント
	Eg[8] = (1 - pa0)*(15*pow(pa3,6) + 45*pow(pa3,4)*pow(pa1,2) + 15*pow(pa3,2)*pow(pa1,4) + pow(pa1,6)) + 15*pa0*pow(pa5,6);						// y_1^6
	Eg[9] = (1 - pa0)/2*((15*pow(pa3,4)*(pa7 + pa1*pa2) + pow(pa1,5)*pa2 + 5*pow(pa1,4)*pa7 + 10*pow(pa1,3)*pa2*pow(pa3,2) + 30*pa7*pow(pa1,2)*pow(pa3,2))
		+ (15*pow(pa3,4)*(pa9 + pa1*pa2) + pow(pa1,5)*pa2 + 5*pow(pa1,4)*pa9 + 10*pow(pa1,3)*pa2*pow(pa3,2) + 30*pa9*pow(pa1,2)*pow(pa3,2)))
		+ 15*pa0*pow(pa5,4)*pa8;									// y_1^5*y_2
	Eg[10] = (1 - pa0)/2*((3*pow(pa3,4)*pow(pa4,2) + 12*pow(pa3*pa7,2)+ 3*pow(pa2,2)*pow(pa3,4) + 6*pow(pa1*pa2*pa3,2) + pow(pa1,4)*pow(pa2,2) + 24*pa1*pa2*pow(pa3,2)*pa7
		+ 8*pow(pa1,3)*pa2*pa7 + 6*pow(pa1*pa3*pa4,2) + 12*pow(pa1*pa7,2) + pow(pa1,4)*pow(pa4,2))
		+ (3*pow(pa3,4)*pow(pa4,2) + 12*pow(pa3*pa9,2)+ 3*pow(pa2,2)*pow(pa3,4) + 6*pow(pa1*pa2*pa3,2) + pow(pa1,4)*pow(pa2,2) + 24*pa1*pa2*pow(pa3,2)*pa9
		+ 8*pow(pa1,3)*pa2*pa9 + 6*pow(pa1*pa3*pa4,2) + 12*pow(pa1*pa9,2) + pow(pa1,4)*pow(pa4,2)))
		+ 3*pa0*(pow(pa5,4)*pow(pa6,2) + 4*pow(pa8,2)*pow(pa5,2));								// y_1^4*y_2^2
	Eg[11] = (1 - pa0)/2*((6*pow((pa7 + pa1*pa2),3) + 9*pa7*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)) + pa1*pa2*(9*pow(pa3*pa4,2)
		+ 3*pow(pa1*pa4,2) + 3*pow(pa2*pa3,2) - 5*pow(pa1*pa2,2)))
		+ (6*pow((pa9 + pa1*pa2),3) + 9*pa9*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)) + pa1*pa2*(9*pow(pa3*pa4,2)
		+ 3*pow(pa1*pa4,2) + 3*pow(pa2*pa3,2) - 5*pow(pa1*pa2,2))))
		+ pa0*(6*pow(pa8,3) + 9*pa8*pow(pa5*pa6,2));															// y_1^3*y_2^3
	Eg[12] = (1 - pa0)/2*((3*(pow(pa4,4)*pow(pa3,2) + 12*pow(pa4*pa7,2) + 3*pow(pa1,2)*pow(pa4,4) + 6*pow(pa1*pa2*pa4,2) + pow(pa2,4)*pow(pa1,2)) + 24*pa1*pa2*pow(pa4,2)*pa7
		+ 8*pow(pa2,3)*pa1*pa7 + 6*pow(pa2*pa3*pa4,2) + 12*pow(pa2*pa7,2)+ pow(pa2,4)*pow(pa3,2))
		+ (3*(pow(pa4,4)*pow(pa3,2) + 12*pow(pa4*pa9,2) + 3*pow(pa1,2)*pow(pa4,4) + 6*pow(pa1*pa2*pa4,2) + pow(pa2,4)*pow(pa1,2)) + 24*pa1*pa2*pow(pa4,2)*pa9
		+ 8*pow(pa2,3)*pa1*pa9 + 6*pow(pa2*pa3*pa4,2) + 12*pow(pa2*pa9,2)+ pow(pa2,4)*pow(pa3,2)))
		+ 3*pa0*(pow(pa6,4)*pow(pa5,2) + 4*pow(pa8,2)*pow(pa6,2));													// y_1^2*y_2^4
	Eg[13] = (1 - pa0)/2*((15*pow(pa4,4)*(pa7 + pa1*pa2) + pow(pa2,5)*pa1 + 5*pow(pa2,4)*pa7 + 10*pow(pa2,3)*pa1*pow(pa4,2) + 30*pa7*pow(pa2,2)*pow(pa4,2))
		+ (15*pow(pa4,4)*(pa9 + pa1*pa2) + pow(pa2,5)*pa1 + 5*pow(pa2,4)*pa9 + 10*pow(pa2,3)*pa1*pow(pa4,2) + 30*pa9*pow(pa2,2)*pow(pa4,2)))
		+ 15*pa0*pow(pa6,4)*pa8;	// y_1*y_2^5
	Eg[14] = (1 - pa0)*(15*pow(pa4,6) + 45*pow(pa4,4)*pow(pa2,2) + 15*pow(pa4,2)*pow(pa2,4) + pow(pa2,6))
		+ 15*pa0*pow(pa6,6);																		// y_2^6
	
	//8次モーメント
	Eg[15] = (1 - pa0)*(pow(pa1,8) + 28*pow(pa3,2)*pow(pa1,6) + 210*pow(pa3*pa1,4) + 420*pow(pa3,6)*pow(pa1,2) + 105*pow(pa3,8))
		+ 105*pa0*pow(pa5,8);																		// y_1^8

	Eg[16] = (1 - pa0)/2*((pow(pa1,7)*pa2 + 21*pow(pa3,2)*pow(pa1,5)*pa2 + 105*pow(pa3,4)*pow(pa1,3)*pa2 + 105*pow(pa3,6)*pa1*pa2 + pa7*(7*pow(pa1,6)
		+ 105*pow(pa3,2)*pow(pa1,4)+ 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6)))
		+ (pow(pa1,7)*pa2 + 21*pow(pa3,2)*pow(pa1,5)*pa2 + 105*pow(pa3,4)*pow(pa1,3)*pa2 + 105*pow(pa3,6)*pa1*pa2 + pa9*(7*pow(pa1,6)
		+ 105*pow(pa3,2)*pow(pa1,4)+ 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6))))
		+ 105*pa0*pa8*pow(pa5,6);																	// y_1^7*y_2
	
	Eg[17] = (1 - pa0)/2*((pow(pa1*pa1*pa1*pa2,2) + 15*pow(pa1*pa1*pa2*pa3,2) + 45*pow(pa1*pa2*pa3*pa3,2) + 15*pow(pa2*pa3*pa3*pa3,2) + 12*pa7*pow(pa1,5)*pa2
		+ 120*pa7*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa7*pa1*pa2*pow(pa3,4) + pow(pa1*pa1*pa1*pa4,2) + 15*pow(pa1*pa1*pa3*pa4,2) + 30*pow(pa7*pa1*pa1,2)
		+ 45*pow(pa1*pa3*pa3*pa4,2) + 180*pow(pa7*pa1*pa3,2) + 15*pow(pa3*pa3*pa3*pa4,2) + 90*pow(pa7*pa3*pa3,2))
		+ (pow(pa1*pa1*pa1*pa2,2) + 15*pow(pa1*pa1*pa2*pa3,2) + 45*pow(pa1*pa2*pa3*pa3,2) + 15*pow(pa2*pa3*pa3*pa3,2) + 12*pa9*pow(pa1,5)*pa2
		+ 120*pa9*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa9*pa1*pa2*pow(pa3,4) + pow(pa1*pa1*pa1*pa4,2) + 15*pow(pa1*pa1*pa3*pa4,2) + 30*pow(pa9*pa1*pa1,2)
		+ 45*pow(pa1*pa3*pa3*pa4,2) + 180*pow(pa9*pa1*pa3,2) + 15*pow(pa3*pa3*pa3*pa4,2) + 90*pow(pa9*pa3*pa3,2)))
		+ 15*pa0*(pow(pa5*pa5*pa5*pa6,2) + 6*pow(pa8*pa5*pa5,2));													// y_1^6y_2^2
	
	Eg[18] = (1 - pa0)/2*(pow(pa1,5)*pow(pa2,3) + 10*pow(pa1*pa2,3)*pow(pa3,2) + 45*pa1*pow(pa2,3)*pow(pa3,4) + 15*pa7*pow(pa1*pa1*pa2,2)
		+ 90*pa7*pow(pa1*pa2*pa3,2) + 45*pa7*pow(pa2*pa3*pa3,2) + 3*pow(pa1,5)*pa2*pow(pa4,2) + 30*pow(pa1,3)*pa2*pow(pa3*pa4,2)
		+ 60*pow(pa7,2)*pow(pa1,3)*pa2 + 45*pa1*pa2*pow(pa3*pa3*pa4,2)+ 180*pow(pa7,2)*pa1*pa2*pow(pa3,2) + 15*pa7*pow(pa1*pa1*pa4,2)
		+ 90*pa7*pow(pa1*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa1,2) + 45*pa7*pow(pa3*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa3,2)
		+ (pow(pa1,5)*pow(pa2,3) + 10*pow(pa1*pa2,3)*pow(pa3,2) + 45*pa1*pow(pa2,3)*pow(pa3,4) + 15*pa9*pow(pa1*pa1*pa2,2)
		+ 90*pa9*pow(pa1*pa2*pa3,2) + 45*pa9*pow(pa2*pa3*pa3,2) + 3*pow(pa1,5)*pa2*pow(pa4,2) + 30*pow(pa1,3)*pa2*pow(pa3*pa4,2)
		+ 60*pow(pa9,2)*pow(pa1,3)*pa2 + 45*pa1*pa2*pow(pa3*pa3*pa4,2)+ 180*pow(pa9,2)*pa1*pa2*pow(pa3,2) + 15*pa9*pow(pa1*pa1*pa4,2)
		+ 90*pa9*pow(pa1*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa1,2) + 45*pa9*pow(pa3*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa3,2)))
		+ pa0*(45*pa8*pow(pa5*pa5*pa6,2) + 60*pow(pa8,3)*pow(pa5,2));													// y_1^5*y_2^3
	
	Eg[19] = (1 - pa0)/2*(pow(pa1*pa2,4) + 6*pow(pa1*pa2*pa2*pa3,2) + 3*pow(pa2*pa3,4) + 16*pa7*pow(pa1*pa2,3) + 48*pa7*pa1*pow(pa2,3)*pow(pa3,2)
		+ 6*pow(pa1*pa1*pa2*pa4,2) + 36*pow(pa1*pa2*pa3*pa4,2) + 72*pow(pa7*pa1*pa2,2) + 18*pow(pa2*pa3*pa3*pa4,2) + 72*pow(pa7*pa2*pa3,2)
		+ 48*pa7*pow(pa1,3)*pa2*pow(pa4,2) + 144*pa7*pa1*pa2*pow(pa3*pa4,2) + 96*pa1*pa2*pow(pa7,3) + 3*pow(pa1*pa4,4) + 18*pow(pa1*pa3*pa4*pa4,2)
		+ 72*pow(pa7*pa1*pa4,2) + 9*pow(pa3*pa4,4) + 72*pow(pa7*pa3*pa4,2) + 24*pow(pa7,4)
		+ (pow(pa1*pa2,4) + 6*pow(pa1*pa2*pa2*pa3,2) + 3*pow(pa2*pa3,4) + 16*pa9*pow(pa1*pa2,3) + 48*pa9*pa1*pow(pa2,3)*pow(pa3,2)
		+ 6*pow(pa1*pa1*pa2*pa4,2) + 36*pow(pa1*pa2*pa3*pa4,2) + 72*pow(pa9*pa1*pa2,2) + 18*pow(pa2*pa3*pa3*pa4,2) + 72*pow(pa9*pa2*pa3,2)
		+ 48*pa9*pow(pa1,3)*pa2*pow(pa4,2) + 144*pa9*pa1*pa2*pow(pa3*pa4,2) + 96*pa1*pa2*pow(pa9,3) + 3*pow(pa1*pa4,4) + 18*pow(pa1*pa3*pa4*pa4,2)
		+ 72*pow(pa9*pa1*pa4,2) + 9*pow(pa3*pa4,4) + 72*pow(pa9*pa3*pa4,2) + 24*pow(pa9,4)))
		+ pa0*(9*pow(pa5*pa6,4) + 72*pow(pa8*pa5*pa6,2) + 24*pow(pa8,4));												// y_1^4*y_2*^4
	
	Eg[20] = (1 - pa0)/2*(pow(pa2,5)*pow(pa1,3) + 10*pow(pa1*pa2,3)*pow(pa4,2) + 45*pa2*pow(pa1,3)*pow(pa4,4) + 15*pa7*pow(pa2*pa2*pa1,2)
		+ 90*pa7*pow(pa1*pa2*pa4,2) + 45*pa7*pow(pa1*pa4*pa4,2) + 3*pow(pa2,5)*pa1*pow(pa3,2) + 30*pow(pa2,3)*pa1*pow(pa3*pa4,2)
		+ 60*pow(pa7,2)*pow(pa2,3)*pa1 + 45*pa1*pa2*pow(pa4*pa4*pa3,2) + 180*pow(pa7,2)*pa1*pa2*pow(pa4,2) + 15*pa7*pow(pa2*pa2*pa3,2)
		+ 90*pa7*pow(pa2*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa2,2) + 45*pa7*pow(pa4*pa4*pa3,2) + 60*pow(pa7,3)*pow(pa4,2)
		+ (pow(pa2,5)*pow(pa1,3) + 10*pow(pa1*pa2,3)*pow(pa4,2) + 45*pa2*pow(pa1,3)*pow(pa4,4) + 15*pa9*pow(pa2*pa2*pa1,2)
		+ 90*pa9*pow(pa1*pa2*pa4,2) + 45*pa9*pow(pa1*pa4*pa4,2) + 3*pow(pa2,5)*pa1*pow(pa3,2) + 30*pow(pa2,3)*pa1*pow(pa3*pa4,2)
		+ 60*pow(pa9,2)*pow(pa2,3)*pa1 + 45*pa1*pa2*pow(pa4*pa4*pa3,2) + 180*pow(pa9,2)*pa1*pa2*pow(pa4,2) + 15*pa9*pow(pa2*pa2*pa3,2)
		+ 90*pa9*pow(pa2*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa2,2) + 45*pa9*pow(pa4*pa4*pa3,2) + 60*pow(pa9,3)*pow(pa4,2)))
		+ pa0*(45*pa8*pow(pa6*pa6*pa5,2) + 60*pow(pa8,3)*pow(pa6,2));													// y_1^3*y_w^5


	double v_result_moment_eq[NUM_OF_MOMENT_EQUATION];
	for (tmp=0; tmp<NUM_OF_MOMENT_EQUATION; tmp++)
	{
		v_result_moment_eq[tmp] = 0.0;
	}

	// モーメント方程式を解く
	gsl_matrix_view m_cf_moment_eq		= gsl_matrix_view_array(cf_moment_eq, NUM_OF_MOMENT_EQUATION, 21);
	gsl_matrix_view m_moment			= gsl_matrix_view_array(Eg, 21, 1);
	gsl_matrix_view m_result_moment_eq	= gsl_matrix_view_array(v_result_moment_eq, NUM_OF_MOMENT_EQUATION, 1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &m_cf_moment_eq.matrix, &m_moment.matrix, 0.0, &m_result_moment_eq.matrix);

	// 計算結果を配列に保存
	double array_result_moment_eq[NUM_OF_MOMENT_EQUATION];
	for (tmp=0; tmp<NUM_OF_MOMENT_EQUATION; tmp++)
	{
		array_result_moment_eq[tmp] = gsl_matrix_get(&m_result_moment_eq.matrix, tmp, 0);	// 先の行列計算の答えを配列にする
	}
	array_result_moment_eq[2]  += dG[1];
	array_result_moment_eq[7]  += dG[3];
	array_result_moment_eq[14] += dG[5];
	
	// 補正係数を含めた結果を格納
	this->_r[NUM_OF_MOMENT_EQUATION]	= {};
	for (tmp=0; tmp<NUM_OF_MOMENT_EQUATION; tmp++)
	{
		this->_r[tmp]	= keisu[tmp]*(array_result_moment_eq[tmp]-y[tmp]);
	}
}

double momentEqFunc::F1()
{
	return this->_r[0];
}

double momentEqFunc::F2()
{
	return this->_r[1];
}

double momentEqFunc::F3()
{
	return this->_r[2];
}

double momentEqFunc::F4()
{
	return this->_r[3];
}

double momentEqFunc::F5()
{
	return this->_r[4];
}

double momentEqFunc::F6()
{
	return this->_r[5];
}

double momentEqFunc::F7()
{
	return this->_r[6];
}

double momentEqFunc::F8()
{
	return this->_r[7];
}

double momentEqFunc::F9()
{
	return this->_r[8];
}

double momentEqFunc::F10()
{
	return this->_r[9];
}

double momentEqFunc::F11()
{
	return this->_r[10];
}

double momentEqFunc::F12()
{
	return this->_r[11];
}

double momentEqFunc::F13()
{
	return this->_r[12];
}

double momentEqFunc::F14()
{
	return this->_r[13];
}

double momentEqFunc::F15()
{
	return this->_r[14];
}