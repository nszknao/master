/**
* メモ
* m_：マトリクス，v_：ベクトル，cf_：係数
**/

struct data
{
	size_t n;	// 方程式の数
	size_t p;	// パラメータの数
	double *y;
	
	//系のパラメータ
	double zeta;
	double epi;
	
	//入力のモーメント4次まで。
	double *dG;
};

int expb_f (const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t NUM_OF_MOMENT_EQUATION	= ((struct data *)data)->n;
	size_t NUM_OF_PARAMETER			= ((struct data *)data)->p;
	double *y						= ((struct data *)data)->y;
	double zeta						= ((struct data *)data)->zeta;
	double epi						= ((struct data *)data)->epi;
	double *dG						= ((struct data *)data)->dG;
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
	keisu[1] = 24./5.;
	keisu[2] = 24./6.;
	keisu[3] = 24./16.;
	keisu[4] = 24./19.;
	keisu[5] = 1.;
	keisu[6] = 24./32.;
	keisu[7] = 24./36.;
	keisu[8] = 24./64.;
	keisu[9] = 24./79.;
	keisu[10] = 24./98.;
	keisu[11] = 24./117.;
	keisu[12] = 24./144;
	keisu[13] = 24./183.;
	keisu[14] = 24./216;

	// カウント変数
	size_t tmp;
	
	double param[NUM_OF_PARAMETER];	// (a, μ1, μ2, σ11, σ12, σ21, σ22, k1, k2, k3)

	for (tmp=0; tmp<NUM_OF_PARAMETER; tmp++)
	{
		param[tmp] = gsl_vector_get(x, tmp);
	}
	// ゆとり作戦
	double pa0 = param[0], pa1 = param[1], pa2 = param[2], pa3 = param[3], pa4 = param[4], pa5 = param[5], pa6 = param[6], pa7 = param[7], pa8 = param[8], pa9 = param[9];

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
	
	// 補正係数を含めた結果をfに格納
	for (tmp=0; tmp<NUM_OF_MOMENT_EQUATION; tmp++)
	{
		gsl_vector_set(f, tmp, keisu[tmp]*(array_result_moment_eq[tmp]-y[tmp]));
	}

	return GSL_SUCCESS;
}

/* ヤコビ行列を定義 */
int expb_df (const gsl_vector * x, void *data, gsl_matrix *J)
{
	size_t NUM_OF_MOMENT_EQUATION	= ((struct data *)data)->n;
	size_t NUM_OF_PARAMETER			= ((struct data *)data)->p;
	double *y						= ((struct data *)data)->y;
	double zeta						= ((struct data *)data)->zeta;
	double epi						= ((struct data *)data)->epi;
	double *dG						= ((struct data *)data)->dG;
	double cf_moment_eq[] =						//モーメント方程式 係数行列 15x21
	{
		0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		-1,-2*zeta,1,-1*epi,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,-2,-4*zeta,0,-2*epi,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	       	
		0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,-1,-2*zeta,3,0,0,-1*epi,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,-2,-4*zeta,2,0,0,-2*epi,0,0,0,0,0,0,0,0,0,0,0,
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
	keisu[1] = 24./5.;
	keisu[2] = 24./6.;
	keisu[3] = 24./16.;
	keisu[4] = 24./19.;
	keisu[5] = 1.;
	keisu[6] = 24./32.;
	keisu[7] = 24./36.;
	keisu[8] = 24./64.;
	keisu[9] = 24./79.;
	keisu[10] = 24./98.;
	keisu[11] = 24./117.;
	keisu[12] = 24./144;
	keisu[13] = 24./183.;
	keisu[14] = 24./216;

	// カウント変数
	size_t tmp, tmp_param, tmp_moment;

	double param[NUM_OF_PARAMETER];
	for(tmp=0; tmp<NUM_OF_PARAMETER; tmp++)
	{
	  param[tmp] = gsl_vector_get(x, tmp);
	}
	// ゆとりパラメータ
	double pa0 = param[0], pa1 = param[1], pa2 = param[2], pa3 = param[3], pa4 = param[4], pa5 = param[5], pa6 = param[6], pa7 = param[7], pa8 = param[8], pa9 = param[9];

	/*************************** ヤコビアン（開始） ***************************/
	double jacoby[21][NUM_OF_PARAMETER];

	// ２次モーメントのヤコビ
	jacoby[0][0] = -(pow(pa3,2) + pow(pa1,2)) + pow(pa5,2);
	jacoby[0][1] = 2*(1 - pa0)*pa1;
	jacoby[0][2] = 0;
	jacoby[0][3] = 2*(1 - pa0)*pa3;
	jacoby[0][4] = 0;
	jacoby[0][5] = 2*pa0*pa5;
	jacoby[0][6] = 0;
	jacoby[0][7] = 0;
	jacoby[0][8] = 0;
	jacoby[0][9] = 0;

	jacoby[1][0] = -1/2*((pa7 + pa1*pa2) + (pa9 + pa1*pa2)) + pa8;
	jacoby[1][1] = (1 - pa0)/2*pa2;
	jacoby[1][2] = (1 - pa0)/2*pa1;
	jacoby[1][3] = 0;
	jacoby[1][4] = 0;
	jacoby[1][5] = 0;
	jacoby[1][6] = 0;
	jacoby[1][7] = (1 - pa0)/2;
	jacoby[1][8] = pa0;
	jacoby[1][9] = (1 - pa0)/2;
				
	jacoby[2][0] = -(pow(pa4,2)+ pow(pa2,2)) + pow(pa6,2);
	jacoby[2][1] = 0;
	jacoby[2][2] = 2*(1 - pa0)*pa2;
	jacoby[2][3] = 0;
	jacoby[2][4] = 2*(1 - pa0)*pa4;
	jacoby[2][5] = 0;
	jacoby[2][6] = 2*pa0*pa6;
	jacoby[2][7] = 0;
	jacoby[2][8] = 0;
	jacoby[2][9] = 0;

	// ４次モーメントのヤコビ
	jacoby[3][0] = -(3*pow(pa3,4) + 6*pow(pa1,2)*pow(pa3,2) + pow(pa1,4)) + 3*pow(pa5,4);
	jacoby[3][1] = (1 - pa0)/2*(12*pa1*pow(pa3,2) + 4*pow(pa1,3));
	jacoby[3][2] = 0;
	jacoby[3][3] = 12*(1 - pa0)*(pow(pa3,3) + pow(pa1,2)*pa3);
	jacoby[3][4] = 0;
	jacoby[3][5] = 12*pa0*pow(pa5,3);
	jacoby[3][6] = 0;
	jacoby[3][7] = 0;
	jacoby[3][8] = 0;
	jacoby[3][9] = 0;
	
	jacoby[4][0] = -1/2*((3*pow(pa3,2)*(pa7 + pa1*pa2) + pow(pa1,2)*(3*pa7 + pa1*pa2)) + (3*pow(pa3,2)*(pa9 + pa1*pa2) + pow(pa1,2)*(3*pa9 + pa1*pa2))) + 3*pow(pa5,2)*pa8;
	jacoby[4][1] = (1 - pa0)/2*((3*pa2*pow(pa3,2) + 6*pa1*pa7 + 3* pow(pa1,2)*pa2) + (3*pa2*pow(pa3,2) + 6*pa1*pa9 + 3* pow(pa1,2)*pa2));
	jacoby[4][2] = (1 - pa0)*(3*pa1*pow(pa3,2) + pow(pa1,3));
	jacoby[4][3] = 3*pa3*(1 - pa0)*((pa7 + pa1*pa2) + (pa9 + pa1*pa2));
	jacoby[4][4] = 0;
	jacoby[4][5] = 6*pa0*pa5*pa8;
	jacoby[4][6] = 0;
	jacoby[4][7] = 3/2*(1 - pa0)*(pow(pa3,2) + pow(pa1,2));
	jacoby[4][8] = 3*pa0*pow(pa5,2);
	jacoby[4][9] = 3/2*(1 - pa0)*(pow(pa3,2) + pow(pa1,2));

	jacoby[5][0] = -1/2*((pow(pa3,2)*pow(pa4,2) + pow(pa1,2)*pow(pa4,2) + pow(pa2,2)*pow(pa3,2) + 2*pow(pa7,2) + 4*pa7*pa1*pa2 + pow(pa1,2)*pow(pa2,2))
			+ (pow(pa3,2)*pow(pa4,2) + pow(pa1,2)*pow(pa4,2) + pow(pa2,2)*pow(pa3,2) + 2*pow(pa9,2) + 4*pa9*pa1*pa2 + pow(pa1,2)*pow(pa2,2)))
			+ pow(pa5,2)*pow(pa6,2) + 2*pow(pa8,2);
	jacoby[5][1] = (1 - pa0)/2*((2*pa1*pow(pa4,2) + 4*pa2*pa7 + 2*pa1*pow(pa2,2)) + (2*pa1*pow(pa4,2) + 4*pa2*pa9 + 2*pa1*pow(pa2,2)));
	jacoby[5][2] = (1 - pa0)/2*((2*pa2*pow(pa3,2) + 4*pa1*pa7 + 2*pa2*pow(pa1,2)) + (2*pa2*pow(pa3,2) + 4*pa1*pa9 + 2*pa2*pow(pa1,2)));
	jacoby[5][3] = (1 - pa0)*(2*pa3*pow(pa4,2) + 2*pow(pa2,2)*pa3);
	jacoby[5][4] = (1 - pa0)*(2*pa4*pow(pa3,2) + 2*pow(pa1,2)*pa4);
	jacoby[5][5] = 2*pa0*pa5*pow(pa6,2);
	jacoby[5][6] = 2*pa0*pa6*pow(pa5,2);
	jacoby[5][7] = 2*(1 - pa0)*(pa7 + pa1*pa2);
	jacoby[5][8] = 4*pa0*pa8;
	jacoby[5][9] = 2*(1 - pa0)*(pa9 + pa1*pa2);

	jacoby[6][0] = -1/2*((3*pow(pa4,2)*(pa7 + pa1*pa2) + pow(pa2,2)*(3*pa7 + pa1*pa2)) + (3*pow(pa4,2)*(pa9 + pa1*pa2) + pow(pa2,2)*(3*pa9 + pa1*pa2))) + 3*pow(pa6,2)*pa8;
	jacoby[6][1] = (1 - pa0)*(3*pa2*pow(pa4,2) + pow(pa2,3));
	jacoby[6][2] = (1 - pa0)/2*((3*pa1*pow(pa4,2) + 6*pa2*pa7 + 3*pa1*pow(pa2,2)) + (3*pa1*pow(pa4,2) + 6*pa2*pa9 + 3*pa1*pow(pa2,2)));
	jacoby[6][3] = 0;
	jacoby[6][4] = 3*pa4*(1 - pa0)*((pa7 + pa1*pa2) + (pa9 + pa1*pa2));
	jacoby[6][5] = 0;
	jacoby[6][6] = 6*pa0*pa6*pa8;
	jacoby[6][7] = 3/2*(1 - pa0)*(pow(pa4,2) + pow(pa2,2));
	jacoby[6][8] = 3*pa0*pow(pa6,2);
	jacoby[6][9] = 3/2*(1 - pa0)*(pow(pa4,2) + pow(pa2,2));

	jacoby[7][0] = -(3*pow(pa4,4) + 6*pow(pa2,2)*pow(pa4,2) + pow(pa2,4)) + 3*pow(pa6,4);
	jacoby[7][1] = 0;
	jacoby[7][2] = (1 - pa0)*(12*pa2*pow(pa4,2) + 4* pow(pa2,3));
	jacoby[7][3] = 0;
	jacoby[7][4] = 12*(1 - pa0)*(pow(pa4,3) + pow(pa2,2)*pa4);
	jacoby[7][5] = 0;
	jacoby[7][6] = 12*pa0*pow(pa6,3);
	jacoby[7][7] = 0;
	jacoby[7][8] = 0;

	// ６次モーメントのヤコビ
	jacoby[8][0] = -(15*pow(pa3,6) + 45*pow(pa3,4)*pow(pa1,2) + 15*pow(pa3,2)*pow(pa1,4) + pow(pa1,6)) + 15*pow(pa5,6);
	jacoby[8][1] = (1 - pa0)*(90*pow(pa3,4)*pa1 + 60*pow(pa3,2)*pow(pa1,3) + 6*pow(pa1,5));
	jacoby[8][2] = 0;
	jacoby[8][3] = (1 - pa0)*(90*pow(pa3,5) + 180*pow(pa1,2)*pow(pa3,3) + 30*pa3*pow(pa1,4));
	jacoby[8][4] = 0;
	jacoby[8][5] = 90*pa0*pow(pa5,5);
	jacoby[8][6] = 0;
	jacoby[8][7] = 0;
	jacoby[8][8] = 0;
	jacoby[8][9] = 0;

	jacoby[9][0] = -1/2*((15*pow(pa3,4)*(pa7 + pa1*pa2) + pow(pa1,5)*pa2 + 5*pow(pa1,4)*pa7 + 10*pow(pa1,3)*pa2*pow(pa3,2) + 30*pa7*pow(pa1,2)*pow(pa3,2))
			+ (15*pow(pa3,4)*(pa9 + pa1*pa2) + pow(pa1,5)*pa2 + 5*pow(pa1,4)*pa9 + 10*pow(pa1,3)*pa2*pow(pa3,2) + 30*pa9*pow(pa1,2)*pow(pa3,2)))
			+ 15*pow(pa5,4)*pa8;
	jacoby[9][1] = (1 - pa0)/2*((15*pa2*pow(pa3,4) + 5*pow(pa1,4)*pa2 + 20*pow(pa1,3)*pa7 + 30*pow(pa1,2)*pa2*pow(pa3,2) + 60*pow(pa3,2)*pa1*pa7)
			+ (15*pa2*pow(pa3,4) + 5*pow(pa1,4)*pa2 + 20*pow(pa1,3)*pa9 + 30*pow(pa1,2)*pa2*pow(pa3,2) + 60*pow(pa3,2)*pa1*pa9));
	jacoby[9][2] = (1 - pa0)*(15*pa1*pow(pa3,4) + pow(pa1,5) + 10*pow(pa1,3)* pow(pa3,2));
	jacoby[9][3] = (1 - pa0)/2*((60*pow(pa3,3)*(pa7 + pa1*pa2) + 20*pow(pa1,3)*pa2*pa3 + 60*pa3*pow(pa1,2)*pa7)
			+ (60*pow(pa3,3)*(pa9 + pa1*pa2) + 20*pow(pa1,3)*pa2*pa3 + 60*pa3*pow(pa1,2)*pa9));
	jacoby[9][4] = 0;
	jacoby[9][5] = 60*pa0*pow(pa5,3)*pa8;
	jacoby[9][6] = 0;
	jacoby[9][7] = (1 - pa0)/2*(15*pow(pa3,4) + 5*pow(pa1,4) + 30*pow(pa1,2)*pow(pa3,2));
	jacoby[9][8] = 15*pa0*pow(pa5,4);
	jacoby[9][9] = (1 - pa0)/2*(15*pow(pa3,4) + 5*pow(pa1,4) + 30*pow(pa1,2)*pow(pa3,2));

	jacoby[10][0] = -1/2*((3*(pow(pa3,4)*pow(pa4,2) + 12*pow(pa3*pa7,2)+ 3*pow(pa2,2)*pow(pa3,4) + 6*pow(pa1*pa2*pa3,2) + pow(pa1,4)*pow(pa2,2))
			+ 24*pa1*pa2*pow(pa3,2)*pa7 + 8*pow(pa1,3)*pa2*pa7 + 6*pow(pa1*pa3*pa4,2) + 12*pow(pa1*pa7,2) + pow(pa1,4)*pow(pa4,2))
			+ (3*(pow(pa3,4)*pow(pa4,2) + 12*pow(pa3*pa9,2)+ 3*pow(pa2,2)*pow(pa3,4) + 6*pow(pa1*pa2*pa3,2) + pow(pa1,4)*pow(pa2,2))
			+ 24*pa1*pa2*pow(pa3,2)*pa9 + 8*pow(pa1,3)*pa2*pa9 + 6*pow(pa1*pa3*pa4,2) + 12*pow(pa1*pa9,2) + pow(pa1,4)*pow(pa4,2)))
			+ 3*(pow(pa5,4)*pow(pa6,2) + 4*pow(pa8*pa5,2));
	jacoby[10][1] = (1 - pa0)/2*((12*pa1*pow(pa2*pa3,2) + 4*pow(pa2,2)*pow(pa1,3) + 24*pow(pa3,2)*pa2*pa7 + 24*pow(pa1,2)*pa2*pa7 + 12*pa1*pow(pa3*pa4,2)
			+ 24*pow(pa7,2)*pa1 + 4*pow(pa1,3)*pow(pa4,2))
			+ (12*pa1*pow(pa2*pa3,2) + 4*pow(pa2,2)*pow(pa1,3) + 24*pow(pa3,2)*pa2*pa9 + 24*pow(pa1,2)*pa2*pa9 + 12*pa1*pow(pa3*pa4,2)
			+ 24*pow(pa9,2)*pa1 + 4*pow(pa1,3)*pow(pa4,2)));
	jacoby[10][2] = (1 - pa0)/2*((12*pow(pa3*pa1,2)*pa2 + 6*pow(pa3,4)*pa2 + 24*pa7*pow(pa3,2)*pa1 + 2*pow(pa1,4)*pa2 + 8*pa7*pow(pa1,3))
			+ (12*pow(pa3*pa1,2)*pa2 + 6*pow(pa3,4)*pa2 + 24*pa9*pow(pa3,2)*pa1 + 2*pow(pa1,4)*pa2 + 8*pa9*pow(pa1,3)));
	jacoby[10][3] = (1 - pa0)/2*((12*pow(pa4,2)*pow(pa3,3) + 24*pa3*pow(pa7,2) + 12*pow(pa2,2)*pow(pa3,3) + 12*pow(pa1*pa2,2)*pa3 + 48*pa1*pa2*pa3*pa7 + 12*pow(pa1*pa4,2)*pa3)
			+ (12*pow(pa4,2)*pow(pa3,3) + 24*pa3*pow(pa9,2) + 12*pow(pa2,2)*pow(pa3,3) + 12*pow(pa1*pa2,2)*pa3 + 48*pa1*pa2*pa3*pa9 + 12*pow(pa1*pa4,2)*pa3));
	jacoby[10][4] = pa4*(1 - pa0)*(6*pow(pa3,4) + 12*pow(pa3*pa1,2) + 2*pow(pa1,4));
	jacoby[10][5] = 3*pa0*(4*pow(pa6,2)*pow(pa5,3) + 8*pow(pa8,2)*pa5);
	jacoby[10][6] = 6*pa0*pow(pa5,4)*pa6;
	jacoby[10][7] = (1 - pa0)/2*(24*(pow(pa3,2)*pa7 + 24*pow(pa3,2)*pa1*pa2 + 24*pow(pa1,2)*pa7 + 8*pow(pa1,3)*pa2));
	jacoby[10][8] = 24*pa0*pow(pa5,2)*pa8;
	jacoby[10][9] = (1 - pa0)/2*(24*(pow(pa3,2)*pa9 + 24*pow(pa3,2)*pa1*pa2 + 24*pow(pa1,2)*pa9 + 8*pow(pa1,3)*pa2));

	jacoby[11][0] = -1/2*((6*pow((pa7+pa1*pa2),3) + 9*pa7*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)) + pa1*pa2*(9*pow(pa3*pa4,2)
			+ 3*pow(pa1*pa4,2) + 3*pow(pa2*pa3,2) - 5*pow(pa1*pa2,2)))
			+ (6*pow((pa9+pa1*pa2),3) + 9*pa9*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)) + pa1*pa2*(9*pow(pa3*pa4,2)
			+ 3*pow(pa1*pa4,2) + 3*pow(pa2*pa3,2) - 5*pow(pa1*pa2,2))))
			+ 6*pow(pa8,3) + 9*pa8*pow(pa5*pa6,2);
	jacoby[11][1] = (1 - pa0)/2*((18*pow((pa7+pa1*pa2),2)*pa2 + 18*pa7*pa1*(pow(pa4,2) - pow(pa2,2)) + pa2*(9*pow(pa3*pa4,2) + 3*pow(pa2*pa3,2))
			+ 9*pa2*pow(pa4*pa1,2) - 15*pow(pa1,2)*pow(pa2,3))
			+ (18*pow((pa9+pa1*pa2),2)*pa2 + 18*pa9*pa1*(pow(pa4,2) - pow(pa2,2)) + pa2*(9*pow(pa3*pa4,2) + 3*pow(pa2*pa3,2))
			+ 9*pa2*pow(pa4*pa1,2) - 15*pow(pa1,2)*pow(pa2,3)));
	jacoby[11][2] = (1 - pa0)/2*((18*pow((pa7+pa1*pa2),2)*pa1 + 18*pa7*pa2*(pow(pa3,2) - pow(pa1,2)) + pa1*(9*pow(pa3*pa4,2) + 3*pow(pa1*pa4,2))
			+ 9*pa1*pow(pa3*pa2,2) - 15*pow(pa2,2)*pow(pa1,3))
			+ (18*pow((pa9+pa1*pa2),2)*pa1 + 18*pa9*pa2*(pow(pa3,2) - pow(pa1,2)) + pa1*(9*pow(pa3*pa4,2) + 3*pow(pa1*pa4,2))
			+ 9*pa1*pow(pa3*pa2,2) - 15*pow(pa2,2)*pow(pa1,3)));
	jacoby[11][3] = (1 - pa0)*pa3*((9*pa7*(pow(pa4,2) + pow(pa2,2)) + pa1*pa2*(9*pow(pa4,2) + 3*pow(pa2,2)))
			+ (9*pa9*(pow(pa4,2) + pow(pa2,2)) + pa1*pa2*(9*pow(pa4,2) + 3*pow(pa2,2))));
	jacoby[11][4] = (1 - pa0)*pa4*((9*pa7*(pow(pa3,2) + pow(pa1,2)) + pa1*pa2*(9*pow(pa3,2) + 3*pow(pa1,2)))
			+ (9*pa9*(pow(pa3,2) + pow(pa1,2)) + pa1*pa2*(9*pow(pa3,2) + 3*pow(pa1,2))));
	jacoby[11][5] = 18*pa0*pa8*pa5*pow(pa6,2);
	jacoby[11][6] = 18*pa0*pa8*pa6*pow(pa5,2);
	jacoby[11][7] = (1 - pa0)/2*(18*pow((pa7+pa1*pa2),2) + 9*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)));
	jacoby[11][8] = pa0*(18*pow(pa8,2) + 9*pow(pa5*pa6,2));
	jacoby[11][9] = (1 - pa0)/2*(18*pow((pa9+pa1*pa2),2) + 9*(pow(pa3*pa4,2) + pow(pa1*pa4,2) + pow(pa2*pa3,2) - pow(pa1*pa2,2)));
	
	jacoby[12][0] = -1/2*((3*(pow(pa4,4)*pow(pa3,2) + 12*pow(pa4*pa7,2)+ 3*pow(pa1,2)*pow(pa4,4) + 6*pow(pa1*pa2*pa4,2) + pow(pa2,4)*pow(pa1,2))
			+ 24*pa1*pa2*pow(pa4,2)*pa7 + 8*pow(pa2,3)*pa1*pa7 + 6*pow(pa2*pa3*pa4,2) + 12*pow(pa2*pa7,2) + pow(pa2,4)*pow(pa3,2))
			+ (3*(pow(pa4,4)*pow(pa3,2) + 12*pow(pa4*pa9,2)+ 3*pow(pa1,2)*pow(pa4,4) + 6*pow(pa1*pa2*pa4,2) + pow(pa2,4)*pow(pa1,2))
			+ 24*pa1*pa2*pow(pa4,2)*pa9 + 8*pow(pa2,3)*pa1*pa9 + 6*pow(pa2*pa3*pa4,2) + 12*pow(pa2*pa9,2) + pow(pa2,4)*pow(pa3,2)))
			+ 3*(pow(pa6,4)*pow(pa5,2) + 4*pow(pa8*pa6,2));
	jacoby[12][1] = (1 - pa0)/2*((12*pow(pa4*pa2,2)*pa1 + 6*pow(pa4,4)*pa1 + 24*pa7*pow(pa4,2)*pa2 + 2*pow(pa2,4)*pa1 + 8*pa7*pow(pa2,3))
			+ (12*pow(pa4*pa2,2)*pa1 + 6*pow(pa4,4)*pa1 + 24*pa9*pow(pa4,2)*pa2 + 2*pow(pa2,4)*pa1 + 8*pa9*pow(pa2,3)));
	jacoby[12][2] = (1 - pa0)/2*((12*pa2*pow(pa1*pa4,2) + 4*pow(pa1,2)*pow(pa2,3) + 24*pow(pa4,2)*pa1*pa7 + 24*pow(pa2,2)*pa1*pa7 + 12*pa2*pow(pa3*pa4,2)
			+ 24*pow(pa7,2)*pa1 + 4*pow(pa2,3)*pow(pa3,2))
			+ (12*pa2*pow(pa1*pa4,2) + 4*pow(pa1,2)*pow(pa2,3) + 24*pow(pa4,2)*pa1*pa9 + 24*pow(pa2,2)*pa1*pa9 + 12*pa2*pow(pa3*pa4,2)
			+ 24*pow(pa9,2)*pa2 + 4*pow(pa2,3)*pow(pa3,2)));
	jacoby[12][3] = pa3*(1 - pa0)*(6*pow(pa4,4) + 12*pow(pa4*pa2,2) + 2*pow(pa2,4));
	jacoby[12][4] = (1 - pa0)/2*((12*pow(pa3,2)*pow(pa4,3) + 24*pa4*pow(pa7,2) + 12*pow(pa1,2)*pow(pa4,3) + 12*pow(pa1*pa2,2)*pa4 + 48*pa1*pa2*pa4*pa7 + 12*pow(pa2*pa3,2)*pa4)
			+ (12*pow(pa3,2)*pow(pa4,3) + 24*pa4*pow(pa9,2) + 12*pow(pa1,2)*pow(pa4,3) + 12*pow(pa1*pa2,2)*pa4 + 48*pa1*pa2*pa4*pa9 + 12*pow(pa2*pa3,2)*pa4));
	jacoby[12][5] = 6*pa0*pow(pa6,4)*pa5;
	jacoby[12][6] = 3*pa0*(4*pow(pa5,2)*pow(pa6,3) + 8*pow(pa8,2)*pa6);
	jacoby[12][7] = (1 - pa0)/2*(24*(pow(pa4,2)*pa7 + 24*pow(pa4,2)*pa1*pa2 + 24*pow(pa2,2)*pa7 + 8*pow(pa2,3)*pa1));
	jacoby[12][8] = 24*pa0*pow(pa6,2)*pa8;
	jacoby[12][9] = (1 - pa0)/2*(24*(pow(pa4,2)*pa9 + 24*pow(pa4,2)*pa1*pa2 + 24*pow(pa2,2)*pa9 + 8*pow(pa2,3)*pa1));

	jacoby[13][0] = -1/2*((15*pow(pa4,4)*(pa7 + pa1*pa2) + pow(pa2,5)*pa1 + 5*pow(pa2,4)*pa7 + 10*pow(pa2,3)*pa1*pow(pa4,2) + 30*pa7*pow(pa2,2)*pow(pa4,2))
			+ (15*pow(pa4,4)*(pa9 + pa1*pa2) + pow(pa2,5)*pa1 + 5*pow(pa2,4)*pa9 + 10*pow(pa2,3)*pa1*pow(pa4,2) + 30*pa9*pow(pa2,2)*pow(pa4,2)))
			+ 15*pow(pa6,4)*pa8;
	jacoby[13][1] = (1 - pa0)*(15*pa2*pow(pa4,4) + pow(pa2,5) + 10*pow(pa2,3)* pow(pa4,2));
	jacoby[13][2] = (1 - pa0)/2*((15*pa1*pow(pa4,4) + 5*pow(pa2,4)*pa1 + 20*pow(pa2,3)*pa7 + 30*pow(pa2,2)*pa1*pow(pa4,2) + 60*pow(pa4,2)*pa2*pa7)
			+ (15*pa1*pow(pa4,4) + 5*pow(pa2,4)*pa1 + 20*pow(pa2,3)*pa9 + 30*pow(pa2,2)*pa1*pow(pa4,2) + 60*pow(pa4,2)*pa2*pa9));
	jacoby[13][3] = 0;
	jacoby[13][4] = (1 - pa0)/2*((60*pow(pa4,3)*(pa7 + pa1*pa2) + 20*pow(pa2,3)*pa1*pa4 + 60*pa4*pow(pa2,2)*pa7)
			+ (60*pow(pa4,3)*(pa9 + pa1*pa2) + 20*pow(pa2,3)*pa1*pa4 + 60*pa4*pow(pa2,2)*pa9));
	jacoby[13][5] = 60*pa0*pow(pa6,3)*pa8;
	jacoby[13][6] = 0;
	jacoby[13][7] = (1 - pa0)/2*(15*pow(pa4,4) + 5*pow(pa2,4) + 30*pow(pa2,2)*pow(pa4,2));
	jacoby[13][8] = 15*pa0*pow(pa6,4);
	jacoby[13][9] = (1 - pa0)/2*(15*pow(pa4,4) + 5*pow(pa2,4) + 30*pow(pa2,2)*pow(pa4,2));
	
	jacoby[14][0] = -(15*pow(pa4,6) + 45*pow(pa4,4)*pow(pa2,2) + 15*pow(pa4,2)*pow(pa2,4) + pow(pa2,6)) + 15*pow(pa6,6);
	jacoby[14][1] = 0;
	jacoby[14][2] = (1 - pa0)*(90*pow(pa4,4)*pa2 + 60*pow(pa4,2)*pow(pa2,3) + 6*pow(pa2,5));
	jacoby[14][3] = 0;
	jacoby[14][4] = (1 - pa0)*(90*pow(pa4,5) + 180*pow(pa2,2)*pow(pa4,3) + 30*pa4*pow(pa2,4));
	jacoby[14][5] = 0;
	jacoby[14][6] = 90*pa0*pow(pa6,5);
	jacoby[14][7] = 0;
	jacoby[14][8] = 0;
	jacoby[14][9] = 0;
	
	// ８次モーメントのヤコビ
	jacoby[15][0] = -(pow(pa1,8) + 28*pow(pa3,2)*pow(pa1,6) + 210*pow(pa3*pa1,4) + 420*pow(pa3,6)*pow(pa1,2) + 105*pow(pa3,8))
			+ 105*pow(pa5,8);
	jacoby[15][1] = (1 - pa0)*(8*pow(pa1,7) + 168*pow(pa3,2)*pow(pa1,5) + 840*pow(pa3,4)*pow(pa1,3) + 840*pow(pa3,6)*pa1);
	jacoby[15][2] = 0;
	jacoby[15][3] = (1 - pa0)*(56*pa3*pow(pa1,6) + 840*pow(pa3,3)*pow(pa1,4) + 2520*pow(pa3,5)*pow(pa1,2) + 840*pow(pa3,7));
	jacoby[15][4] = 0;
	jacoby[15][5] = 840*pa0*pow(pa5,7);
	jacoby[15][6] = 0;
	jacoby[15][7] = 0;
	jacoby[15][8] = 0;
	
	jacoby[16][0] = -1/2*((pow(pa1,7)*pa2 + 21*pow(pa3,2)*pow(pa1,5)*pa2 + 105*pow(pa3,4)*pow(pa1,3)*pa2 + 105*pow(pa3,6)*pa1*pa2 + pa7*(7*pow(pa1,6)
			+ 105*pow(pa3,2)*pow(pa1,4) + 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6)))
			+ (pow(pa1,7)*pa2 + 21*pow(pa3,2)*pow(pa1,5)*pa2 + 105*pow(pa3,4)*pow(pa1,3)*pa2 + 105*pow(pa3,6)*pa1*pa2 + pa9*(7*pow(pa1,6)
			+ 105*pow(pa3,2)*pow(pa1,4) + 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6))))
			+ 105*pa0*pa8*pow(pa5,6);
	jacoby[16][1] = (1 - pa0)/2*((7*pow(pa1,6)*pa2 + 105*pow(pa3,2)*pow(pa1,4)*pa2 + 315*pow(pa3,4)*pow(pa1,2)*pa2 + 105*pow(pa3,6)*pa2 + pa7*(42*pow(pa1,5)
			+ 420*pow(pa3,2)*pow(pa1,3) + 630*pow(pa3,4)*pa1))
			+ (7*pow(pa1,6)*pa2 + 105*pow(pa3,2)*pow(pa1,4)*pa2 + 315*pow(pa3,4)*pow(pa1,2)*pa2 + 105*pow(pa3,6)*pa2 + pa9*(42*pow(pa1,5)
			+ 420*pow(pa3,2)*pow(pa1,3) + 630*pow(pa3,4)*pa1)));
	jacoby[16][2] = (1 - pa0)*(pow(pa1,7) + 21*pow(pa3,2)*pow(pa1,5) + 105*pow(pa3,4)*pow(pa1,3) + 105*pow(pa3,6)*pa1);
	jacoby[16][3] = (1 - pa0)/2*((42*pa3*pow(pa1,5)*pa2 + 420*pow(pa3,3)*pow(pa1,3)*pa2 + 630*pow(pa3,5)*pa1*pa2 + pa7*(210*pa3*pow(pa1,4) + 1260*pow(pa3,3)*pow(pa1,2) + 630*pow(pa3,5)))
			+ (42*pa3*pow(pa1,5)*pa2 + 420*pow(pa3,3)*pow(pa1,3)*pa2 + 630*pow(pa3,5)*pa1*pa2 + pa9*(210*pa3*pow(pa1,4) + 1260*pow(pa3,3)*pow(pa1,2) + 630*pow(pa3,5))));
	jacoby[16][4] = 0;
	jacoby[16][5] = 630*pa0*pa8*pow(pa5,5);
	jacoby[16][6] = 0;
	jacoby[16][7] = (1 - pa0)/2*(7*pow(pa1,6) + 105*pow(pa3,2)*pow(pa1,4) + 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6));
	jacoby[16][8] = 105*pa0*pow(pa5,6);
	jacoby[16][9] = (1 - pa0)/2*(7*pow(pa1,6) + 105*pow(pa3,2)*pow(pa1,4) + 315*pow(pa3,4)*pow(pa1,2) + 105*pow(pa3,6));

	jacoby[17][0] = -1/2*((pow(pa1*pa1*pa1*pa2,2) + 15*pow(pa1*pa1*pa2*pa3,2) + 45*pow(pa1*pa2*pa3*pa3,2) + 15*pow(pa2*pa3*pa3*pa3,2) + 12*pa7*pow(pa1,5)*pa2
			+ 120*pa7*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa7*pa1*pa2*pow(pa3,4) + pow(pa1*pa1*pa1*pa4,2) + 15*pow(pa1*pa1*pa3*pa4,2)
			+ 30*pow(pa7*pa1*pa1,2) + 45*pow(pa1*pa3*pa3*pa4,2) + 180*pow(pa7*pa1*pa3,2) + 15*pow(pa3*pa3*pa3*pa4,2)+ 90*pow(pa7*pa3*pa3,2))
			+ (pow(pa1*pa1*pa1*pa2,2) + 15*pow(pa1*pa1*pa2*pa3,2) + 45*pow(pa1*pa2*pa3*pa3,2) + 15*pow(pa2*pa3*pa3*pa3,2) + 12*pa9*pow(pa1,5)*pa2
			+ 120*pa9*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa9*pa1*pa2*pow(pa3,4) + pow(pa1*pa1*pa1*pa4,2) + 15*pow(pa1*pa1*pa3*pa4,2)
			+ 30*pow(pa9*pa1*pa1,2) + 45*pow(pa1*pa3*pa3*pa4,2) + 180*pow(pa9*pa1*pa3,2) + 15*pow(pa3*pa3*pa3*pa4,2) + 90*pow(pa9*pa3*pa3,2)))
			+ 15*(pow(pa5*pa5*pa5*pa6,2) + 6*pow(pa8*pa5*pa5,2));
	jacoby[17][1] = (1 - pa0)/2*((6*pow(pa1,5)*pow(pa2,2) + 60*pow(pa1,3)*pow(pa2*pa3,2) + 90*pa1*pow(pa2*pa3*pa3,2)+ 60*pa7*pow(pa1,4)*pa2 + 360*pa7*pow(pa1*pa3,2)*pa2
			+ 180*pa7*pa2*pow(pa3,4) + 6*pow(pa1,5)*pow(pa4,2)+ 60*pow(pa3*pa4,2)*pow(pa1,3) + 120*pow(pa7,2)*pow(pa1,3) + 90*pa1*pow(pa3*pa3*pa4,2) + 360*pow(pa3*pa7,2)*pa1)
			+ (6*pow(pa1,5)*pow(pa2,2) + 60*pow(pa1,3)*pow(pa2*pa3,2) + 90*pa1*pow(pa2*pa3*pa3,2)+ 60*pa9*pow(pa1,4)*pa2 + 360*pa9*pow(pa1*pa3,2)*pa2
			+ 180*pa9*pa2*pow(pa3,4) + 6*pow(pa1,5)*pow(pa4,2)+ 60*pow(pa3*pa4,2)*pow(pa1,3) + 120*pow(pa9,2)*pow(pa1,3) + 90*pa1*pow(pa3*pa3*pa4,2) + 360*pow(pa3*pa9,2)*pa1));
	jacoby[17][2] = (1 - pa0)/2*((2*pow(pa1,6)*pa2 + 30*pow(pa1*pa1*pa3,2)*pa2 + 90*pa2*pow(pa1*pa3*pa3,2) + 30*pa2*pow(pa3,6) + 12*pa7*pow(pa1,5) + 120*pa7*pow(pa1,3)*pow(pa3,2) + 180*pa7*pa1*pow(pa3,4))
			+ (2*pow(pa1,6)*pa2 + 30*pow(pa1*pa1*pa3,2)*pa2 + 90*pa2*pow(pa1*pa3*pa3,2) + 30*pa2*pow(pa3,6) + 12*pa9*pow(pa1,5) + 120*pa9*pow(pa1,3)*pow(pa3,2) + 180*pa9*pa1*pow(pa3,4)));
	jacoby[17][3] = (1 - pa0)/2*((30*pow(pa1*pa1*pa2,2)*pa3 + 180*pow(pa1*pa2,2)*pow(pa3,3) + 90*pow(pa2,2)*pow(pa3,5)+ 240*pa7*pow(pa1,3)*pa2*pa3 + 720*pa7*pa1*pa2*pow(pa3,3)
			+ 30*pow(pa1*pa1*pa4,2)*pa3 + 180*pow(pa1*pa4,2)*pow(pa3,3) + 360*pow(pa7*pa1,2)*pa3 + 90*pow(pa4,2)*pow(pa3,5) + 360*pow(pa7,2)*pow(pa3,3))
			+ (30*pow(pa1*pa1*pa2,2)*pa3 + 180*pow(pa1*pa2,2)*pow(pa3,3) + 90*pow(pa2,2)*pow(pa3,5) + 240*pa9*pow(pa1,3)*pa2*pa3 + 720*pa9*pa1*pa2*pow(pa3,3)
			+ 30*pow(pa1*pa1*pa4,2)*pa3 + 180*pow(pa1*pa4,2)*pow(pa3,3)+ 360*pow(pa9*pa1,2)*pa3 + 90*pow(pa4,2)*pow(pa3,5) + 360*pow(pa9,2)*pow(pa3,3)));
	jacoby[17][4] = (1 - pa0)*(2*pow(pa1,6)*pa4 + 30*pow(pa1*pa1*pa3,2)*pa4 + 90*pow(pa1*pa3*pa3,2)*pa4 + 30*pow(pa3,6)*pa4);
	jacoby[17][5] = pa0*(90*pow(pa5,5)*pow(pa6,2) + 360*pow(pa8,2)*pow(pa6,3));
	jacoby[17][6] = pa0*30*pow(pa5,6)*pa6;
	jacoby[17][7] = (1 - pa0)/2*(12*pow(pa1,5)*pa2 + 120*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa1*pa2*pow(pa3,4) + 60*pa7*pow(pa1,4) + 360*pa7*pow(pa1*pa3,2) + 180*pa7*pow(pa3,4));
	jacoby[17][8] = pa0*180*pa8*pow(pa5,4);
	jacoby[17][9] = (1 - pa0)/2*(12*pow(pa1,5)*pa2 + 120*pow(pa1,3)*pa2*pow(pa3,2) + 180*pa1*pa2*pow(pa3,4) + 60*pa9*pow(pa1,4) + 360*pa9*pow(pa1*pa3,2) + 180*pa9*pow(pa3,4));
	
	jacoby[18][0] = -1/2*((pow(pa1,5)*pow(pa2,3) + 10*pow(pa1*pa2,3)*pow(pa3,2) + 45*pa1*pow(pa2,3)*pow(pa3,4) + 15*pa7*pow(pa1*pa1*pa2,2)
			+ 90*pa7*pow(pa1*pa2*pa3,2) + 45*pa7*pow(pa2*pa3*pa3,2) + 3*pow(pa1,5)*pa2*pow(pa4,2) + 30*pow(pa1,3)*pa2*pow(pa3*pa4,2)
			+ 60*pow(pa7,2)*pow(pa1,3)*pa2 + 45*pa1*pa2*pow(pa3*pa3*pa4,2) + 180*pow(pa7,2)*pa1*pa2*pow(pa3,2) + 15*pa7*pow(pa1*pa1*pa4,2)
			+ 90*pa7*pow(pa1*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa1,2) + 45*pa7*pow(pa3*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa3,2))
			+ (pow(pa1,5)*pow(pa2,3) + 10*pow(pa1*pa2,3)*pow(pa3,2) + 45*pa1*pow(pa2,3)*pow(pa3,4) + 15*pa9*pow(pa1*pa1*pa2,2)
			+ 90*pa9*pow(pa1*pa2*pa3,2) + 45*pa9*pow(pa2*pa3*pa3,2) + 3*pow(pa1,5)*pa2*pow(pa4,2) + 30*pow(pa1,3)*pa2*pow(pa3*pa4,2)
			+ 60*pow(pa9,2)*pow(pa1,3)*pa2 + 45*pa1*pa2*pow(pa3*pa3*pa4,2) + 180*pow(pa9,2)*pa1*pa2*pow(pa3,2) + 15*pa9*pow(pa1*pa1*pa4,2)
			+ 90*pa9*pow(pa1*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa1,2) + 45*pa9*pow(pa3*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa3,2)))
			+ (45*pa8*pow(pa5*pa5*pa6,2) + 60*pow(pa8,3)*pow(pa5,2));
	jacoby[18][1] = (1 - pa0)/2*((5*pow(pa1,4)*pow(pa2,3) + 30*pow(pa1*pa3,2)*pow(pa2,3) + 45*pow(pa2,3)*pow(pa3,4) + 60*pa7*pow(pa1,3)*pow(pa2,2) + 180*pa7*pa1*pow(pa2*pa3,2)
			+ 15*pow(pa1*pa1*pa4,2)*pa2 + 90*pow(pa1*pa3*pa4,2)*pa2+ 180*pow(pa7*pa1,2)*pa2 + 45*pa2*pow(pa3*pa3*pa4,2) + 180*pow(pa7*pa3,2)*pa2 + 60*pa7*pow(pa1,3)*pow(pa4,2)
			+ 180*pa7*pa1*pow(pa3*pa4,2) + 120*pow(pa7,3)*pa1)
			+ (5*pow(pa1,4)*pow(pa2,3) + 30*pow(pa1*pa3,2)*pow(pa2,3) + 45*pow(pa2,3)*pow(pa3,4) + 60*pa9*pow(pa1,3)*pow(pa2,2) + 180*pa9*pa1*pow(pa2*pa3,2)
			+ 15*pow(pa1*pa1*pa4,2)*pa2 + 90*pow(pa1*pa3*pa4,2)*pa2+ 180*pow(pa9*pa1,2)*pa2 + 45*pa2*pow(pa3*pa3*pa4,2) + 180*pow(pa9*pa3,2)*pa2 + 60*pa9*pow(pa1,3)*pow(pa4,2)
			+ 180*pa9*pa1*pow(pa3*pa4,2) + 120*pow(pa9,3)*pa1));
	jacoby[18][2] = (1 - pa0)/2*((3*pow(pa1,5)*pow(pa2,2) + 30*pow(pa1,3)*pow(pa2*pa3,2) + 135*pa1*pow(pa2*pa3*pa3,2) + 30*pa7*pow(pa1,4)*pa2 + 180*pa2*pa7*pow(pa1*pa3,2)
			+ 90*pa7*pa2*pow(pa3,4) + 3*pow(pa1,5)*pow(pa4,2)+ 30*pow(pa1,3)*pow(pa3*pa4,2) + 60*pow(pa7,2)*pow(pa1,3) + 45*pa1*pow(pa3*pa3*pa4,2) + 180*pow(pa7*pa3,2)*pa1)
			+ (3*pow(pa1,5)*pow(pa2,2) + 30*pow(pa1,3)*pow(pa2*pa3,2) + 135*pa1*pow(pa2*pa3*pa3,2) + 30*pa9*pow(pa1,4)*pa2 + 180*pa2*pa9*pow(pa1*pa3,2)
			+ 90*pa9*pa2*pow(pa3,4) + 3*pow(pa1,5)*pow(pa4,2)+ 30*pow(pa1,3)*pow(pa3*pa4,2) + 60*pow(pa9,2)*pow(pa1,3) + 45*pa1*pow(pa3*pa3*pa4,2) + 180*pow(pa9*pa3,2)*pa1));
	jacoby[18][3] = (1 - pa0)/2*((20*pow(pa1*pa2,3)*pa3 + 180*pa1*pow(pa2*pa3,3) + 180*pa7*pow(pa1*pa2,2)*pa3 + 180*pa7*pow(pa2,2)*pow(pa3,3) + 60*pow(pa1,3)*pa2*pow(pa4,2)*pa3
			+ 180*pa1*pa2*pow(pa4,2)*pow(pa3,3) + 360*pow(pa7,2)*pa1*pa2*pa3 + 180*pa7*pow(pa1*pa4,2)*pa3 + 180*pa7*pow(pa3,3)*pow(pa4,2) + 120*pow(pa7,3)*pa3)
			+ (20*pow(pa1*pa2,3)*pa3 + 180*pa1*pow(pa2*pa3,3) + 180*pa9*pow(pa1*pa2,2)*pa3 + 180*pa9*pow(pa2,2)*pow(pa3,3) + 60*pow(pa1,3)*pa2*pow(pa4,2)*pa3
			+ 180*pa1*pa2*pow(pa4,2)*pow(pa3,3) + 360*pow(pa9,2)*pa1*pa2*pa3 + 180*pa9*pow(pa1*pa4,2)*pa3 + 180*pa9*pow(pa3,3)*pow(pa4,2) + 120*pow(pa9,3)*pa3));
	jacoby[18][4] = (1 - pa0)/2*((6*pow(pa1,5)*pa2*pa4 + 60*pow(pa1,3)*pa2*pow(pa3,2)*pa4 + 90*pa1*pa2*pow(pa3,4)*pa4+ 30*pa7*pow(pa1,4)*pa4 + 180*pa7*pow(pa1*pa3,2)*pa4 + 90*pa7*pow(pa3,4)*pa4)
			+ (6*pow(pa1,5)*pa2*pa4 + 60*pow(pa1,3)*pa2*pow(pa3,2)*pa4 + 90*pa1*pa2*pow(pa3,4)*pa4+ 30*pa9*pow(pa1,4)*pa4 + 180*pa9*pow(pa1*pa3,2)*pa4 + 90*pa9*pow(pa3,4)*pa4));
	jacoby[18][5] = pa0*(180*pa8*pow(pa5,3)*pow(pa6,2) + 120*pow(pa8,3)*pa5);
	jacoby[18][6] = pa0*90*pa8*pow(pa5,4)*pa6;
	jacoby[18][7] = (1 - pa0)/2*(15*pow(pa1*pa1*pa2,2) + 90*pow(pa1*pa2*pa3,2) + 45*pow(pa2,2)*pow(pa3,4) + 120*pa7*pow(pa1,3)*pa2 + 360*pa7*pa1*pa2*pow(pa3,2)
			+ 15*pow(pa1*pa1*pa4,2) + 90*pow(pa1*pa3*pa4,2) + 180*pow(pa7*pa1,2) + 45*pow(pa3,4)*pow(pa4,2) + 180*pow(pa7*pa3,2));
	jacoby[18][8] = pa0*(45*pow(pa5*pa5*pa6,2) + 180*pow(pa8*pa5,2));
	jacoby[18][9] = (1 - pa0)/2*(15*pow(pa1*pa1*pa2,2) + 90*pow(pa1*pa2*pa3,2) + 45*pow(pa2,2)*pow(pa3,4) + 120*pa9*pow(pa1,3)*pa2 + 360*pa9*pa1*pa2*pow(pa3,2)
			+ 15*pow(pa1*pa1*pa4,2) + 90*pow(pa1*pa3*pa4,2) + 180*pow(pa9*pa1,2) + 45*pow(pa3,4)*pow(pa4,2) + 180*pow(pa9*pa3,2));

	jacoby[19][0] = -1/2* ((pow(pa1*pa2,4) + 6*pow(pa1*pa2*pa2*pa3,2) + 3*pow(pa2*pa3,4) + 16*pa7*pow(pa1*pa2,3) + 48*pa7*pa1*pow(pa2,3)*pow(pa3,2) + 6*pow(pa1*pa1*pa2*pa4,2)
			+ 36*pow(pa1*pa2*pa3*pa4,2) + 72*pow(pa7*pa1*pa2,2) + 18*pow(pa2*pa3*pa3*pa4,2) + 72*pow(pa7*pa2*pa3,2)+ 48*pa7*pow(pa1,3)*pa2*pow(pa4,2) + 144*pa7*pa1*pa2*pow(pa3*pa4,2)
			+ 96*pa1*pa2*pow(pa7,3) + 3*pow(pa1*pa4,4) + 18*pow(pa1*pa3*pa4*pa4,2) + 72*pow(pa7*pa1*pa4,2) + 9*pow(pa3*pa4,4) + 72*pow(pa7*pa3*pa4,2) + 24*pow(pa7,4))
			+ (pow(pa1*pa2,4) + 6*pow(pa1*pa2*pa2*pa3,2)+ 3*pow(pa2*pa3,4) + 16*pa9*pow(pa1*pa2,3) + 48*pa9*pa1*pow(pa2,3)*pow(pa3,2)+ 6*pow(pa1*pa1*pa2*pa4,2)
			+ 36*pow(pa1*pa2*pa3*pa4,2) + 72*pow(pa9*pa1*pa2,2) + 18*pow(pa2*pa3*pa3*pa4,2) + 72*pow(pa9*pa2*pa3,2)+ 48*pa9*pow(pa1,3)*pa2*pow(pa4,2) + 144*pa9*pa1*pa2*pow(pa3*pa4,2)
			+ 96*pa1*pa2*pow(pa9,3) + 3*pow(pa1*pa4,4)+ 18*pow(pa1*pa3*pa4*pa4,2) + 72*pow(pa9*pa1*pa4,2) + 9*pow(pa3*pa4,4) + 72*pow(pa9*pa3*pa4,2) + 24*pow(pa9,4)))
			+ 9*pow(pa5*pa6,4) + 72*pow(pa8*pa5*pa6,2) + 24*pow(pa8,4);
	jacoby[19][1] = (1 - pa0)*((4*pow(pa1,3)*pow(pa2,4) + 12*pa1*pow(pa2*pa2*pa3,2) + 48*pa7*pow(pa2,3)*pow(pa1,2) + 48*pa7*pow(pa2,3)*pow(pa3,2) + 24*pow(pa1,3)*pow(pa2*pa4,2)
			+ 72*pa1*pow(pa2*pa3*pa4,2) + 144*pow(pa7,2)*pow(pa2,2)*pa1+ 144*pa7*pow(pa1*pa4,2)*pa2 + 144*pa7*pa2*pow(pa3*pa4,2) + 96*pa2*pow(pa7,3) + 12*pow(pa4,4)*pow(pa1,3)
			+ 36*pow(pa3*pa4*pa4,2)*pa1 + 144*pow(pa7*pa4,2)*pa1)
			+ (4*pow(pa1,3)*pow(pa2,4) + 12*pa1*pow(pa2*pa2*pa3,2) + 48*pa9*pow(pa2,3)*pow(pa1,2) + 48*pa9*pow(pa2,3)*pow(pa3,2) + 24*pow(pa1,3)*pow(pa2*pa4,2)
			+ 72*pa1*pow(pa2*pa3*pa4,2) + 144*pow(pa9,2)*pow(pa2,2)*pa1+ 144*pa9*pow(pa1*pa4,2)*pa2 + 144*pa9*pa2*pow(pa3*pa4,2) + 96*pa2*pow(pa9,3) + 12*pow(pa4,4)*pow(pa1,3)
			+ 36*pow(pa3*pa4*pa4,2)*pa1 + 144*pow(pa9*pa4,2)*pa1));
	jacoby[19][2] = (1 - pa0)*((4*pow(pa1,4)*pow(pa2,3) + 24*pow(pa1*pa3,2)*pow(pa2,3) + 12*pow(pa2,3)*pow(pa3,4) + 48*pa7*pow(pa1,3)*pow(pa2,2) + 144*pa7*pa1*pow(pa2*pa3,2)
			+ 12*pow(pa1*pa1*pa4,2)*pa2 + 72*pow(pa1*pa3*pa4,2)*pa2 + 144*pow(pa7*pa1,2)*pa2 + 36*pow(pa3*pa3*pa4,2)*pa2 + 144*pow(pa7*pa3,2)*pa2 + 48*pa7*pow(pa1,3)*pow(pa4,2)
			+ 144*pa7*pa1*pow(pa3*pa4,2) + 96*pa1*pow(pa7,3))
			+ (4*pow(pa1,4)*pow(pa2,3) + 24*pow(pa1*pa3,2)*pow(pa2,3) + 12*pow(pa2,3)*pow(pa3,4) + 48*pa9*pow(pa1,3)*pow(pa2,2) + 144*pa9*pa1*pow(pa2*pa3,2)
			+ 12*pow(pa1*pa1*pa4,2)*pa2 + 72*pow(pa1*pa3*pa4,2)*pa2 + 144*pow(pa9*pa1,2)*pa2 + 36*pow(pa3*pa3*pa4,2)*pa2 + 144*pow(pa9*pa3,2)*pa2 + 48*pa9*pow(pa1,3)*pow(pa4,2)
			+ 144*pa9*pa1*pow(pa3*pa4,2) + 96*pa1*pow(pa9,3)));
	jacoby[19][3] = (1 - pa0)*((12*pow(pa1*pa2*pa2,2)*pa3 + 12*pow(pa2,4)*pow(pa3,3) + 96*pa7*pa1*pow(pa2,3)*pa3 + 72*pow(pa1*pa2*pa4,2)*pa3 + 72*pow(pa2*pa4,2)*pow(pa3,3)
			+ 144*pow(pa7*pa2,2)*pa3 + 288*pa7*pa1*pa2*pa3*pow(pa4,2) + 36*pow(pa1*pa4*pa4,2)*pa3+ 36*pow(pa3,3)*pow(pa4,4) + 144*pow(pa7*pa4,2)*pa3)
			+ (12*pow(pa1*pa2*pa2,2)*pa3 + 12*pow(pa2,4)*pow(pa3,3) + 96*pa9*pa1*pow(pa2,3)*pa3 + 72*pow(pa1*pa2*pa4,2)*pa3 + 72*pow(pa2*pa4,2)*pow(pa3,3)
			+ 144*pow(pa9*pa2,2)*pa3 + 288*pa9*pa1*pa2*pa3*pow(pa4,2) + 36*pow(pa1*pa4*pa4,2)*pa3+ 36*pow(pa3,3)*pow(pa4,4) + 144*pow(pa9*pa4,2)*pa3));
	jacoby[19][4] = (1 - pa0)*((12*pow(pa1*pa1*pa2,2)*pa4 + 72*pow(pa1*pa2*pa3,2)*pa4 + 36*pow(pa2*pa3*pa3,2)*pa4 + 96*pa7*pow(pa1,3)*pa2*pa4 + 288*pa7*pa1*pa2*pow(pa3,2)*pa4
			+ 12*pow(pa1,4)*pow(pa4,3) + 72*pow(pa1*pa3,2)*pow(pa4,3) + 144*pow(pa7*pa1,2)*pa4 + 36*pow(pa3,4)*pow(pa4,3) + 144*pow(pa7*pa3,2)*pa4)
			+ (12*pow(pa1*pa1*pa2,2)*pa4 + 72*pow(pa1*pa2*pa3,2)*pa4 + 36*pow(pa2*pa3*pa3,2)*pa4 + 96*pa9*pow(pa1,3)*pa2*pa4 + 288*pa9*pa1*pa2*pow(pa3,2)*pa4
			+ 12*pow(pa1,4)*pow(pa4,3) + 72*pow(pa1*pa3,2)*pow(pa4,3) + 144*pow(pa9*pa1,2)*pa4 + 36*pow(pa3,4)*pow(pa4,3) + 144*pow(pa9*pa3,2)*pa4));
	jacoby[19][5] = pa0*(36*pow(pa5,3)*pow(pa6,4) + 144*pow(pa8*pa6,2)*pa5);
	jacoby[19][6] = pa0*(36*pow(pa5,4)*pow(pa6,3) + 144*pow(pa7*pa5,2)*pa6);
	jacoby[19][7] = (1 - pa0)/2*(16*pow(pa1*pa2,3) + 48*pa1*pow(pa2,3)*pow(pa3,2) + 144*pa7*pow(pa1*pa2,2) + 144*pa7*pow(pa2*pa3,2) + 48*pow(pa1,3)*pa2*pow(pa4,2)
			+ 144*pa1*pa2*pow(pa3*pa4,2) + 144*pa7*pow(pa1*pa4,2) + 288*pa1*pa2*pow(pa7,2) + 144*pa7*pow(pa3*pa4,2) + 96*pow(pa7,3));
	jacoby[19][8] = pa0*(144*pa8*pow(pa5*pa6,2) + 96*pow(pa8,3));
	jacoby[19][9] = (1 - pa0)/2*(16*pow(pa1*pa2,3) + 48*pa1*pow(pa2,3)*pow(pa3,2) + 144*pa9*pow(pa1*pa2,2) + 144*pa9*pow(pa2*pa3,2) + 48*pow(pa1,3)*pa2*pow(pa4,2)
			+ 144*pa1*pa2*pow(pa3*pa4,2) + 144*pa9*pow(pa1*pa4,2) + 288*pa1*pa2*pow(pa9,2) + 144*pa9*pow(pa3*pa4,2) + 96*pow(pa9,3));
	
	jacoby[20][0] = -1/2*((pow(pa2,5)*pow(pa1,3) + 10*pow(pa2*pa1,3)*pow(pa4,2) + 45*pa2*pow(pa1,3)*pow(pa4,4) + 15*pa7*pow(pa2*pa2*pa1,2)
			+ 90*pa7*pow(pa1*pa2*pa4,2) + 45*pa7*pow(pa3*pa2*pa2,2) + 3*pow(pa2,5)*pa1*pow(pa3,2) + 30*pow(pa2,3)*pa1*pow(pa3*pa4,2)
			+ 60*pow(pa7,2)*pow(pa2,3)*pa1 + 45*pa1*pa2*pow(pa4*pa4*pa3,2) + 180*pow(pa7,2)*pa1*pa2*pow(pa4,2) + 15*pa7*pow(pa2*pa2*pa3,2)
			+ 90*pa7*pow(pa2*pa3*pa4,2) + 60*pow(pa7,3)*pow(pa2,2) + 45*pa7*pow(pa4*pa4*pa3,2) + 60*pow(pa7,3)*pow(pa4,2))
			+ (pow(pa2,5)*pow(pa1,3) + 10*pow(pa1*pa2,3)*pow(pa4,2) + 45*pa2*pow(pa1,3)*pow(pa4,4) + 15*pa9*pow(pa2*pa2*pa1,2)
			+ 90*pa9*pow(pa1*pa2*pa4,2) + 45*pa9*pow(pa1*pa4*pa4,2) + 3*pow(pa2,5)*pa1*pow(pa3,2) + 30*pow(pa2,3)*pa1*pow(pa3*pa4,2)
			+ 60*pow(pa9,2)*pow(pa2,3)*pa1 + 45*pa1*pa2*pow(pa4*pa3*pa4,2) + 180*pow(pa9,2)*pa1*pa2*pow(pa4,2) + 15*pa9*pow(pa2*pa2*pa3,2)
			+ 90*pa9*pow(pa2*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa2,2) + 45*pa9*pow(pa4*pa3*pa4,2) + 60*pow(pa9,3)*pow(pa4,2)))
			+ (45*pa8*pow(pa6*pa5*pa6,2) + 60*pow(pa8,3)*pow(pa6,2));
	jacoby[20][1] = (1 - pa0)/2*((5*pow(pa2,4)*pow(pa1,3) + 30*pow(pa2*pa4,2)*pow(pa1,3) + 45*pow(pa1,3)*pow(pa4,4) + 60*pa7*pow(pa2,3)*pow(pa1,2) + 180*pa7*pa2*pow(pa1*pa4,2)
			+ 15*pow(pa2*pa2*pa3,2)*pa1 + 90*pow(pa2*pa3*pa4,2)*pa1+ 180*pow(pa7*pa2,2)*pa1 + 45*pa1*pow(pa4*pa3*pa4,2) + 180*pow(pa7*pa4,2)*pa1 + 60*pa7*pow(pa2,3)*pow(pa3,2)
			+ 180*pa7*pa2*pow(pa3*pa4,2) + 120*pow(pa7,3)*pa2)
			+ (5*pow(pa2,4)*pow(pa1,3) + 30*pow(pa2*pa4,2)*pow(pa1,3) + 45*pow(pa1,3)*pow(pa4,4) + 60*pa9*pow(pa2,3)*pow(pa1,2) + 180*pa9*pa2*pow(pa1*pa4,2)
			+ 15*pow(pa2*pa2*pa3,2)*pa1 + 90*pow(pa2*pa3*pa4,2)*pa1+ 180*pow(pa9*pa2,2)*pa1 + 45*pa1*pow(pa4*pa3*pa4,2) + 180*pow(pa9*pa4,2)*pa1 + 60*pa9*pow(pa2,3)*pow(pa3,2)
			+ 180*pa9*pa2*pow(pa3*pa4,2) + 120*pow(pa9,3)*pa2));
	jacoby[20][2] = (1 - pa0)/2*((3*pow(pa2,5)*pow(pa1,2) + 30*pow(pa2,3)*pow(pa1*pa4,2) + 135*pa2*pow(pa1*pa4*pa4,2) + 30*pa7*pow(pa2,4)*pa1 + 180*pa1*pa7*pow(pa2*pa4,2)
			+ 90*pa7*pa1*pow(pa4,4) + 3*pow(pa2,5)*pow(pa3,2)+ 30*pow(pa2,3)*pow(pa3*pa4,2) + 60*pow(pa7,2)*pow(pa2,3) + 45*pa2*pow(pa4*pa3*pa4,2) + 180*pow(pa7*pa4,2)*pa2)
			+ (3*pow(pa2,5)*pow(pa1,2) + 30*pow(pa2,3)*pow(pa1*pa4,2) + 135*pa2*pow(pa1*pa4*pa4,2) + 30*pa9*pow(pa2,4)*pa1 + 180*pa1*pa9*pow(pa2*pa4,2)
			+ 90*pa9*pa1*pow(pa4,4) + 3*pow(pa2,5)*pow(pa3,2)+ 30*pow(pa2,3)*pow(pa3*pa4,2) + 60*pow(pa9,2)*pow(pa2,3) + 45*pa2*pow(pa4*pa3*pa4,2) + 180*pow(pa9*pa4,2)*pa2));
	jacoby[20][3] = (1 - pa0)/2*((20*pow(pa1*pa2,3)*pa4 + 180*pa2*pow(pa1*pa4,3) + 180*pa7*pow(pa1*pa2,2)*pa4 + 180*pa7*pow(pa1,2)*pow(pa4,3) + 60*pow(pa2,3)*pa1*pow(pa3,2)*pa4
			+ 180*pa1*pa2*pow(pa3,2)*pow(pa4,3) + 360*pow(pa7,2)*pa1*pa2*pa4 + 180*pa7*pow(pa2*pa3,2)*pa4 + 180*pa7*pow(pa4,3)*pow(pa3,2) + 120*pow(pa7,3)*pa4)
			+ (20*pow(pa1*pa2,3)*pa4 + 180*pa2*pow(pa1*pa4,3) + 180*pa9*pow(pa1*pa2,2)*pa4 + 180*pa9*pow(pa1,2)*pow(pa4,3) + 60*pow(pa2,3)*pa1*pow(pa3,2)*pa4
			+ 180*pa1*pa2*pow(pa3,2)*pow(pa4,3) + 360*pow(pa9,2)*pa1*pa2*pa4 + 180*pa9*pow(pa2*pa3,2)*pa4 + 180*pa9*pow(pa4,3)*pow(pa3,2) + 120*pow(pa9,3)*pa4));
	jacoby[20][4] = (1 - pa0)/2*((6*pow(pa2,5)*pa1*pa3 + 60*pow(pa2,3)*pa1*pow(pa4,2)*pa3 + 90*pa1*pa2*pow(pa4,4)*pa3+ 30*pa7*pow(pa2,4)*pa3 + 180*pa7*pow(pa2*pa4,2)*pa3 + 90*pa7*pow(pa4,4)*pa3)
			+ (6*pow(pa2,5)*pa1*pa3 + 60*pow(pa2,3)*pa1*pow(pa4,2)*pa3 + 90*pa1*pa2*pow(pa4,4)*pa3+ 30*pa9*pow(pa2,4)*pa3 + 180*pa9*pow(pa2*pa4,2)*pa4 + 90*pa9*pow(pa4,4)*pa3));
	jacoby[20][5] = pa0*(180*pa8*pow(pa6,3)*pow(pa5,2) + 120*pow(pa8,3)*pa6);
	jacoby[20][6] = pa0*90*pa8*pow(pa6,4)*pa5;
	jacoby[20][7] = (1 - pa0)/2*(15*pow(pa2*pa1*pa2,2) + 90*pow(pa1*pa2*pa4,2) + 45*pow(pa1,2)*pow(pa4,4) + 120*pa7*pow(pa2,3)*pa1 + 360*pa7*pa1*pa2*pow(pa4,2)
			+ 15*pow(pa2*pa2*pa3,2) + 90*pow(pa2*pa3*pa4,2) + 180*pow(pa7*pa2,2) + 45*pow(pa4,4)*pow(pa3,2) + 180*pow(pa7*pa4,2));
	jacoby[20][8] = pa0*(45*pow(pa6*pa5*pa6,2) + 180*pow(pa8*pa6,2));
	jacoby[20][9] = (1 - pa0)/2*(15*pow(pa2*pa1*pa2,2) + 90*pow(pa1*pa2*pa4,2) + 45*pow(pa1,2)*pow(pa4,4) + 120*pa9*pow(pa2,3)*pa1 + 360*pa9*pa1*pa2*pow(pa4,2)
			+ 15*pow(pa2*pa2*pa3,2) + 90*pow(pa2*pa3*pa4,2) + 180*pow(pa9*pa2,2) + 45*pow(pa4,4)*pow(pa3,2) + 180*pow(pa9*pa4,2));

	
	
	// モーメントのヤコビアンをベクトルに変換
	double v_jacoby_moment[210];
	for (tmp_moment=0; tmp_moment<21; tmp_moment++)
	{
		for (tmp_param=0; tmp_param<NUM_OF_PARAMETER; tmp_param++)
		{
			v_jacoby_moment[NUM_OF_PARAMETER*tmp_moment + tmp_param] = jacoby[tmp_moment][tmp_param] ;
		}
	}

	double cf_jacoby_moment_eq[] =
	{
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,

		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
	};

	// モーメント方程式のヤコビアンを計算
	gsl_matrix_view m_cf_jacoby_moment_eq	= gsl_matrix_view_array(cf_jacoby_moment_eq, NUM_OF_MOMENT_EQUATION, NUM_OF_PARAMETER);
	gsl_matrix_view m_cf_moment_eq			= gsl_matrix_view_array(cf_moment_eq, NUM_OF_MOMENT_EQUATION, 21);
	gsl_matrix_view m_jacoby_moment			= gsl_matrix_view_array(v_jacoby_moment, 21, NUM_OF_PARAMETER);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &m_cf_moment_eq.matrix, &m_jacoby_moment.matrix, 0.0, &m_cf_jacoby_moment_eq.matrix);

	// 計算結果をJに格納
	for(tmp_moment=0; tmp_moment<NUM_OF_MOMENT_EQUATION; tmp_moment++)
	{
		for(tmp_param=0; tmp_param<10; tmp_param++)
		{
			gsl_matrix_set(J, tmp_moment, tmp_param, keisu[tmp_moment]*gsl_matrix_get(&m_cf_jacoby_moment_eq.matrix,tmp_moment,tmp_param));
		}
	}

	return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	expb_f(x, data, f);
	expb_df(x, data, J);
	
	return GSL_SUCCESS;
}
