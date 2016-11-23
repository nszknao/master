/**
 * メモ
 * m_：マトリクス，v_：ベクトル，cf_：係数
 */
#include "../include/expfit.h"
#include "../include/common.h"

/**
 * @fn GSLの非線形最小二乗法で使うモーメント方程式の関数値ベクトル
 * @param vector<double> x 方程式のパラメータ (a, μ1, μ2, σ11, σ12, σ21, σ22, k1, k2, k3)
 * @param vector<double> dF 系や入力のパラメータ
 * @param vector<std::size_t> &f 必要な目的関数の番号
 */
std::vector<double>
MomentEq::expb_f (const std::vector<double> &x, const std::vector<double> &dG, const std::vector<std::size_t> &f)
{
    std::size_t i;

	std::vector<double> cf_moment_eq =						// モーメント方程式 係数行列 15x21
	{
		0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		-1,-2*Common::ZETA,1,-1*Common::EPSILON,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,-2,-4*Common::ZETA,0,-2*Common::EPSILON,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	       	
		0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,-1,-2*Common::ZETA,3,0,0,-1*Common::EPSILON,0,0,0,0,0,0,0,0,0,0,0,0,
		dG[1],0,0,0,-2,-4*Common::ZETA,2,0,0,-2*Common::EPSILON,0,0,0,0,0,0,0,0,0,0,0, 
		0,3*dG[1],0,0,0,-3,-6*Common::ZETA,1,0,0,-3*Common::EPSILON,0,0,0,0,0,0,0,0,0,0,
		0,0,6*dG[1],0,0,0,-4,-8*Common::ZETA,0,0,0,-4*Common::EPSILON,0,0,0,0,0,0,0,0,0,

		0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,-1,-2*Common::ZETA,5,0,0,0,0,-1*Common::EPSILON,0,0,0,0,0,
		0,0,0,dG[1],0,0,0,0,0,-2,-4*Common::ZETA,4,0,0,0,0,-2*Common::EPSILON,0,0,0,0,
		0,0,0,0,3*dG[1],0,0,0,0,0,-3,-6*Common::ZETA,3,0,0,0,0,-3*Common::EPSILON,0,0,0,
		dG[3],0,0,0,0,6*dG[1],0,0,0,0,0,-4,-8*Common::ZETA,2,0,0,0,0,-4*Common::EPSILON,0,0,
		0,5*dG[3],0,0,0,0,10*dG[1],0,0,0,0,0,-5,-10*Common::ZETA,1,0,0,0,0,-5*Common::EPSILON,0,
		0,0,15*dG[3],0,0,0,0,15*dG[1],0,0,0,0,0,-6,-12*Common::ZETA,0,0,0,0,0,-6*Common::EPSILON
	};

	// モーメント方程式を計算
	std::vector<double> Eg(Common::NUM_OF_MOMENT);
	MomentEq::getMomentFromParameter(x, Eg);
	// モーメント法手式の結果を保存するvector
	std::vector<double> v_result_moment_eq(Common::NUM_OF_MOMENTEQ, 0.);
	// モーメント方程式を解く
	gsl_matrix_view m_cf_moment_eq		= gsl_matrix_view_array(&cf_moment_eq[0], Common::NUM_OF_MOMENTEQ, Common::NUM_OF_MOMENT);
	gsl_matrix_view m_moment			= gsl_matrix_view_array(&Eg[0], Common::NUM_OF_MOMENT, 1);
	gsl_matrix_view m_result_moment_eq	= gsl_matrix_view_array(&v_result_moment_eq[0], Common::NUM_OF_MOMENTEQ, 1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &m_cf_moment_eq.matrix, &m_moment.matrix, 0.0, &m_result_moment_eq.matrix);
	// 計算結果を配列に保存
	std::vector<double> array_result_moment_eq(Common::NUM_OF_MOMENTEQ);
	for (i=0; i<Common::NUM_OF_MOMENTEQ; ++i) {
		array_result_moment_eq[i] = gsl_matrix_get(&m_result_moment_eq.matrix, i, 0);	// 先の行列計算の答えを配列にする
	}
	array_result_moment_eq[2]  += dG[1];
	array_result_moment_eq[7]  += dG[3];
	array_result_moment_eq[14] += dG[5];
 	
    // 補正係数
	std::vector<double> k(Common::NUM_OF_MOMENTEQ);
	k[0] = 1.26;	k[1] = 0.66;	k[2] = 1.35;	k[3] = 3.81;	k[4] = 2.53;	k[5] = 2.85;	k[6] = 5.16;	k[7] = 7.25;
	k[8] = 19.6;	k[9] = 11.5;	k[10] = 13.7;	k[11] = 18.3;	k[12] = 26.9;	k[13] = 27.9;	k[14] = 318.;
	// 指定した目的関数の値をセット
    std::vector<double> obj;
    obj.reserve(f.size());
    for (i = 0; i < Common::NUM_OF_MOMENTEQ; ++i) {
        if (i == f[i]) {
            //obj.push_back(1./k[i]array_result_moment_eq[i]);
            obj.push_back(array_result_moment_eq[i]);
        }
    }

    return obj;
}

/**
 * @fn パラメータ値からモーメント値を得る
 * @param const std::vector<double> &p パラメータ値
 * @param const std::vector<double> &m モーメント値（領域確保済み）
 */
void MomentEq::getMomentFromParameter(const std::vector<double> &p, std::vector<double> &m)
{
	// ２次モーメント
	m[0] = (1 - p[0])*(pow(p[3],2) + pow(p[1],2)) + p[0]*pow(p[5],2);	// y_1^2
	// m[1] = (1 - p[0])/2*((p[7] + p[1]*p[2]) + (p[9] + p[1]*p[2])) + p[0]*p[8];	// y_1*y_2
	m[1] = 0.0000001;	// y_1*y_2
	m[2] = (1 - p[0])*(pow(p[4],2) + pow(p[2],2)) + p[0]*pow(p[6],2);	// y_2^2

	// ４次モーメント
	m[3] = (1 - p[0])*(3*pow(p[3],4) + 6*pow(p[1],2)*pow(p[3],2) + pow(p[1],4)) + 3*p[0]*pow(p[5],4);	// y_1^4
	// m[4] = (1 - p[0])/2*((3*pow(p[3],2)*(p[7] + p[1]*p[2]) + pow(p[1],2)*(3*p[7] + p[1]*p[2]))
	// 	+ (3*pow(p[3],2)*(p[9] + p[1]*p[2]) + pow(p[1],2)*(3*p[9] + p[1]*p[2]))) + 3*p[0]*pow(p[5],2)*p[8];	// y_1^3*y_2
	m[4] = 0.0000001;	// y_1^3*y_2	
	m[5] = (1 - p[0])/2*((pow(p[3],2)*pow(p[4],2) + pow(p[1],2)*pow(p[4],2) + pow(p[2],2)*pow(p[3],2) + 2*pow(p[7],2) + 4*p[7]*p[1]*p[2] + pow(p[1],2)*pow(p[2],2))
		+ (pow(p[3],2)*pow(p[4],2) + pow(p[1],2)*pow(p[4],2) + pow(p[2],2)*pow(p[3],2) + 2*pow(p[9],2) + 4*p[9]*p[1]*p[2] + pow(p[1],2)*pow(p[2],2)))
		+ p[0]*(pow(p[5],2)*pow(p[6],2) + 2*pow(p[8],2));	// y_1^2*y_2^2
	m[6] = (1 - p[0])/2*((3*pow(p[4],2)*(p[7] + p[1]*p[2]) + pow(p[2],2)*(3*p[7] + p[1]*p[2]))
		+ (3*pow(p[4],2)*(p[9] + p[1]*p[2]) + pow(p[2],2)*(3*p[9] + p[1]*p[2])))
		+ 3*p[0]*pow(p[6],2)*p[8];	// y_1*y_2^3
	m[7] = (1 - p[0])*(3*pow(p[4],4) + 6*pow(p[2],2)*pow(p[4],2) + pow(p[2],4)) + 3*p[0]*pow(p[6],4);	// y_2^4
	
	// ６次モーメント
	m[8] = (1 - p[0])*(15*pow(p[3],6) + 45*pow(p[3],4)*pow(p[1],2) + 15*pow(p[3],2)*pow(p[1],4) + pow(p[1],6)) + 15*p[0]*pow(p[5],6);						// y_1^6
	// m[9] = (1 - p[0])/2*((15*pow(p[3],4)*(p[7] + p[1]*p[2]) + pow(p[1],5)*p[2] + 5*pow(p[1],4)*p[7] + 10*pow(p[1],3)*p[2]*pow(p[3],2) + 30*p[7]*pow(p[1],2)*pow(p[3],2))
	// 	+ (15*pow(p[3],4)*(p[9] + p[1]*p[2]) + pow(p[1],5)*p[2] + 5*pow(p[1],4)*p[9] + 10*pow(p[1],3)*p[2]*pow(p[3],2) + 30*p[9]*pow(p[1],2)*pow(p[3],2)))
	// 	+ 15*p[0]*pow(p[5],4)*p[8];	// y_1^5*y_2
	m[9] = 0.0000001;	// y_1^5*y_2
	m[10] = (1 - p[0])/2*((3*pow(p[3],4)*pow(p[4],2) + 12*pow(p[3]*p[7],2)+ 3*pow(p[2],2)*pow(p[3],4) + 6*pow(p[1]*p[2]*p[3],2) + pow(p[1],4)*pow(p[2],2) + 24*p[1]*p[2]*pow(p[3],2)*p[7]
		+ 8*pow(p[1],3)*p[2]*p[7] + 6*pow(p[1]*p[3]*p[4],2) + 12*pow(p[1]*p[7],2) + pow(p[1],4)*pow(p[4],2))
		+ (3*pow(p[3],4)*pow(p[4],2) + 12*pow(p[3]*p[9],2)+ 3*pow(p[2],2)*pow(p[3],4) + 6*pow(p[1]*p[2]*p[3],2) + pow(p[1],4)*pow(p[2],2) + 24*p[1]*p[2]*pow(p[3],2)*p[9]
		+ 8*pow(p[1],3)*p[2]*p[9] + 6*pow(p[1]*p[3]*p[4],2) + 12*pow(p[1]*p[9],2) + pow(p[1],4)*pow(p[4],2)))
		+ 3*p[0]*(pow(p[5],4)*pow(p[6],2) + 4*pow(p[8],2)*pow(p[5],2));	// y_1^4*y_2^2
	m[11] = (1 - p[0])/2*((6*pow((p[7] + p[1]*p[2]),3) + 9*p[7]*(pow(p[3]*p[4],2) + pow(p[1]*p[4],2) + pow(p[2]*p[3],2) - pow(p[1]*p[2],2)) + p[1]*p[2]*(9*pow(p[3]*p[4],2)
		+ 3*pow(p[1]*p[4],2) + 3*pow(p[2]*p[3],2) - 5*pow(p[1]*p[2],2)))
		+ (6*pow((p[9] + p[1]*p[2]),3) + 9*p[9]*(pow(p[3]*p[4],2) + pow(p[1]*p[4],2) + pow(p[2]*p[3],2) - pow(p[1]*p[2],2)) + p[1]*p[2]*(9*pow(p[3]*p[4],2)
		+ 3*pow(p[1]*p[4],2) + 3*pow(p[2]*p[3],2) - 5*pow(p[1]*p[2],2))))
		+ p[0]*(6*pow(p[8],3) + 9*p[8]*pow(p[5]*p[6],2));	// y_1^3*y_2^3
	m[12] = (1 - p[0])/2*((3*(pow(p[4],4)*pow(p[3],2) + 12*pow(p[4]*p[7],2) + 3*pow(p[1],2)*pow(p[4],4) + 6*pow(p[1]*p[2]*p[4],2) + pow(p[2],4)*pow(p[1],2)) + 24*p[1]*p[2]*pow(p[4],2)*p[7]
		+ 8*pow(p[2],3)*p[1]*p[7] + 6*pow(p[2]*p[3]*p[4],2) + 12*pow(p[2]*p[7],2)+ pow(p[2],4)*pow(p[3],2))
		+ (3*(pow(p[4],4)*pow(p[3],2) + 12*pow(p[4]*p[9],2) + 3*pow(p[1],2)*pow(p[4],4) + 6*pow(p[1]*p[2]*p[4],2) + pow(p[2],4)*pow(p[1],2)) + 24*p[1]*p[2]*pow(p[4],2)*p[9]
		+ 8*pow(p[2],3)*p[1]*p[9] + 6*pow(p[2]*p[3]*p[4],2) + 12*pow(p[2]*p[9],2)+ pow(p[2],4)*pow(p[3],2)))
		+ 3*p[0]*(pow(p[6],4)*pow(p[5],2) + 4*pow(p[8],2)*pow(p[6],2));	// y_1^2*y_2^4
	m[13] = (1 - p[0])/2*((15*pow(p[4],4)*(p[7] + p[1]*p[2]) + pow(p[2],5)*p[1] + 5*pow(p[2],4)*p[7] + 10*pow(p[2],3)*p[1]*pow(p[4],2) + 30*p[7]*pow(p[2],2)*pow(p[4],2))
		+ (15*pow(p[4],4)*(p[9] + p[1]*p[2]) + pow(p[2],5)*p[1] + 5*pow(p[2],4)*p[9] + 10*pow(p[2],3)*p[1]*pow(p[4],2) + 30*p[9]*pow(p[2],2)*pow(p[4],2)))
		+ 15*p[0]*pow(p[6],4)*p[8];	// y_1*y_2^5
	m[14] = (1 - p[0])*(15*pow(p[4],6) + 45*pow(p[4],4)*pow(p[2],2) + 15*pow(p[4],2)*pow(p[2],4) + pow(p[2],6))
		+ 15*p[0]*pow(p[6],6);	// y_2^6
	
	//8次モーメント
	m[15] = (1 - p[0])*(pow(p[1],8) + 28*pow(p[3],2)*pow(p[1],6) + 210*pow(p[3]*p[1],4) + 420*pow(p[3],6)*pow(p[1],2) + 105*pow(p[3],8))
		+ 105*p[0]*pow(p[5],8);	// y_1^8

	m[16] = (1 - p[0])/2*((pow(p[1],7)*p[2] + 21*pow(p[3],2)*pow(p[1],5)*p[2] + 105*pow(p[3],4)*pow(p[1],3)*p[2] + 105*pow(p[3],6)*p[1]*p[2] + p[7]*(7*pow(p[1],6)
		+ 105*pow(p[3],2)*pow(p[1],4)+ 315*pow(p[3],4)*pow(p[1],2) + 105*pow(p[3],6)))
		+ (pow(p[1],7)*p[2] + 21*pow(p[3],2)*pow(p[1],5)*p[2] + 105*pow(p[3],4)*pow(p[1],3)*p[2] + 105*pow(p[3],6)*p[1]*p[2] + p[9]*(7*pow(p[1],6)
		+ 105*pow(p[3],2)*pow(p[1],4)+ 315*pow(p[3],4)*pow(p[1],2) + 105*pow(p[3],6))))
		+ 105*p[0]*p[8]*pow(p[5],6);	// y_1^7*y_2
	
	m[17] = (1 - p[0])/2*((pow(p[1]*p[1]*p[1]*p[2],2) + 15*pow(p[1]*p[1]*p[2]*p[3],2) + 45*pow(p[1]*p[2]*p[3]*p[3],2) + 15*pow(p[2]*p[3]*p[3]*p[3],2) + 12*p[7]*pow(p[1],5)*p[2]
		+ 120*p[7]*pow(p[1],3)*p[2]*pow(p[3],2) + 180*p[7]*p[1]*p[2]*pow(p[3],4) + pow(p[1]*p[1]*p[1]*p[4],2) + 15*pow(p[1]*p[1]*p[3]*p[4],2) + 30*pow(p[7]*p[1]*p[1],2)
		+ 45*pow(p[1]*p[3]*p[3]*p[4],2) + 180*pow(p[7]*p[1]*p[3],2) + 15*pow(p[3]*p[3]*p[3]*p[4],2) + 90*pow(p[7]*p[3]*p[3],2))
		+ (pow(p[1]*p[1]*p[1]*p[2],2) + 15*pow(p[1]*p[1]*p[2]*p[3],2) + 45*pow(p[1]*p[2]*p[3]*p[3],2) + 15*pow(p[2]*p[3]*p[3]*p[3],2) + 12*p[9]*pow(p[1],5)*p[2]
		+ 120*p[9]*pow(p[1],3)*p[2]*pow(p[3],2) + 180*p[9]*p[1]*p[2]*pow(p[3],4) + pow(p[1]*p[1]*p[1]*p[4],2) + 15*pow(p[1]*p[1]*p[3]*p[4],2) + 30*pow(p[9]*p[1]*p[1],2)
		+ 45*pow(p[1]*p[3]*p[3]*p[4],2) + 180*pow(p[9]*p[1]*p[3],2) + 15*pow(p[3]*p[3]*p[3]*p[4],2) + 90*pow(p[9]*p[3]*p[3],2)))
		+ 15*p[0]*(pow(p[5]*p[5]*p[5]*p[6],2) + 6*pow(p[8]*p[5]*p[5],2));	// y_1^6y_2^2
	
	m[18] = (1 - p[0])/2*(pow(p[1],5)*pow(p[2],3) + 10*pow(p[1]*p[2],3)*pow(p[3],2) + 45*p[1]*pow(p[2],3)*pow(p[3],4) + 15*p[7]*pow(p[1]*p[1]*p[2],2)
		+ 90*p[7]*pow(p[1]*p[2]*p[3],2) + 45*p[7]*pow(p[2]*p[3]*p[3],2) + 3*pow(p[1],5)*p[2]*pow(p[4],2) + 30*pow(p[1],3)*p[2]*pow(p[3]*p[4],2)
		+ 60*pow(p[7],2)*pow(p[1],3)*p[2] + 45*p[1]*p[2]*pow(p[3]*p[3]*p[4],2)+ 180*pow(p[7],2)*p[1]*p[2]*pow(p[3],2) + 15*p[7]*pow(p[1]*p[1]*p[4],2)
		+ 90*p[7]*pow(p[1]*p[3]*p[4],2) + 60*pow(p[7],3)*pow(p[1],2) + 45*p[7]*pow(p[3]*p[3]*p[4],2) + 60*pow(p[7],3)*pow(p[3],2)
		+ (pow(p[1],5)*pow(p[2],3) + 10*pow(p[1]*p[2],3)*pow(p[3],2) + 45*p[1]*pow(p[2],3)*pow(p[3],4) + 15*p[9]*pow(p[1]*p[1]*p[2],2)
		+ 90*p[9]*pow(p[1]*p[2]*p[3],2) + 45*p[9]*pow(p[2]*p[3]*p[3],2) + 3*pow(p[1],5)*p[2]*pow(p[4],2) + 30*pow(p[1],3)*p[2]*pow(p[3]*p[4],2)
		+ 60*pow(p[9],2)*pow(p[1],3)*p[2] + 45*p[1]*p[2]*pow(p[3]*p[3]*p[4],2)+ 180*pow(p[9],2)*p[1]*p[2]*pow(p[3],2) + 15*p[9]*pow(p[1]*p[1]*p[4],2)
		+ 90*p[9]*pow(p[1]*p[3]*p[4],2) + 60*pow(p[9],3)*pow(p[1],2) + 45*p[9]*pow(p[3]*p[3]*p[4],2) + 60*pow(p[9],3)*pow(p[3],2)))
		+ p[0]*(45*p[8]*pow(p[5]*p[5]*p[6],2) + 60*pow(p[8],3)*pow(p[5],2));	// y_1^5*y_2^3
	
	m[19] = (1 - p[0])/2*(pow(p[1]*p[2],4) + 6*pow(p[1]*p[2]*p[2]*p[3],2) + 3*pow(p[2]*p[3],4) + 16*p[7]*pow(p[1]*p[2],3) + 48*p[7]*p[1]*pow(p[2],3)*pow(p[3],2)
		+ 6*pow(p[1]*p[1]*p[2]*p[4],2) + 36*pow(p[1]*p[2]*p[3]*p[4],2) + 72*pow(p[7]*p[1]*p[2],2) + 18*pow(p[2]*p[3]*p[3]*p[4],2) + 72*pow(p[7]*p[2]*p[3],2)
		+ 48*p[7]*pow(p[1],3)*p[2]*pow(p[4],2) + 144*p[7]*p[1]*p[2]*pow(p[3]*p[4],2) + 96*p[1]*p[2]*pow(p[7],3) + 3*pow(p[1]*p[4],4) + 18*pow(p[1]*p[3]*p[4]*p[4],2)
		+ 72*pow(p[7]*p[1]*p[4],2) + 9*pow(p[3]*p[4],4) + 72*pow(p[7]*p[3]*p[4],2) + 24*pow(p[7],4)
		+ (pow(p[1]*p[2],4) + 6*pow(p[1]*p[2]*p[2]*p[3],2) + 3*pow(p[2]*p[3],4) + 16*p[9]*pow(p[1]*p[2],3) + 48*p[9]*p[1]*pow(p[2],3)*pow(p[3],2)
		+ 6*pow(p[1]*p[1]*p[2]*p[4],2) + 36*pow(p[1]*p[2]*p[3]*p[4],2) + 72*pow(p[9]*p[1]*p[2],2) + 18*pow(p[2]*p[3]*p[3]*p[4],2) + 72*pow(p[9]*p[2]*p[3],2)
		+ 48*p[9]*pow(p[1],3)*p[2]*pow(p[4],2) + 144*p[9]*p[1]*p[2]*pow(p[3]*p[4],2) + 96*p[1]*p[2]*pow(p[9],3) + 3*pow(p[1]*p[4],4) + 18*pow(p[1]*p[3]*p[4]*p[4],2)
		+ 72*pow(p[9]*p[1]*p[4],2) + 9*pow(p[3]*p[4],4) + 72*pow(p[9]*p[3]*p[4],2) + 24*pow(p[9],4)))
		+ p[0]*(9*pow(p[5]*p[6],4) + 72*pow(p[8]*p[5]*p[6],2) + 24*pow(p[8],4));	// y_1^4*y_2*^4
	
	m[20] = (1 - p[0])/2*(pow(p[2],5)*pow(p[1],3) + 10*pow(p[1]*p[2],3)*pow(p[4],2) + 45*p[2]*pow(p[1],3)*pow(p[4],4) + 15*p[7]*pow(p[2]*p[2]*p[1],2)
		+ 90*p[7]*pow(p[1]*p[2]*p[4],2) + 45*p[7]*pow(p[1]*p[4]*p[4],2) + 3*pow(p[2],5)*p[1]*pow(p[3],2) + 30*pow(p[2],3)*p[1]*pow(p[3]*p[4],2)
		+ 60*pow(p[7],2)*pow(p[2],3)*p[1] + 45*p[1]*p[2]*pow(p[4]*p[4]*p[3],2) + 180*pow(p[7],2)*p[1]*p[2]*pow(p[4],2) + 15*p[7]*pow(p[2]*p[2]*p[3],2)
		+ 90*p[7]*pow(p[2]*p[3]*p[4],2) + 60*pow(p[7],3)*pow(p[2],2) + 45*p[7]*pow(p[4]*p[4]*p[3],2) + 60*pow(p[7],3)*pow(p[4],2)
		+ (pow(p[2],5)*pow(p[1],3) + 10*pow(p[1]*p[2],3)*pow(p[4],2) + 45*p[2]*pow(p[1],3)*pow(p[4],4) + 15*p[9]*pow(p[2]*p[2]*p[1],2)
		+ 90*p[9]*pow(p[1]*p[2]*p[4],2) + 45*p[9]*pow(p[1]*p[4]*p[4],2) + 3*pow(p[2],5)*p[1]*pow(p[3],2) + 30*pow(p[2],3)*p[1]*pow(p[3]*p[4],2)
		+ 60*pow(p[9],2)*pow(p[2],3)*p[1] + 45*p[1]*p[2]*pow(p[4]*p[4]*p[3],2) + 180*pow(p[9],2)*p[1]*p[2]*pow(p[4],2) + 15*p[9]*pow(p[2]*p[2]*p[3],2)
		+ 90*p[9]*pow(p[2]*p[3]*p[4],2) + 60*pow(p[9],3)*pow(p[2],2) + 45*p[9]*pow(p[4]*p[4]*p[3],2) + 60*pow(p[9],3)*pow(p[4],2)))
		+ p[0]*(45*p[8]*pow(p[6]*p[6]*p[5],2) + 60*pow(p[8],3)*pow(p[6],2));	// y_1^3*y_w^5
}
