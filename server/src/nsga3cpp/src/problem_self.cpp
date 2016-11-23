#include "problem_self.h"
#include "alg_individual.h"
#include "../../../include/expfit.h"

using namespace std;

CProblemSelf::CProblemSelf(std::size_t num_vars, std::size_t num_objs, double labmda, double alpha):
	BProblem("Self"),
	num_vars_(num_vars),
	num_objs_(num_objs),
	lambda_(labmda),
	alpha_(alpha)
{
	// define the domain of variables here
	lbs_.resize(num_vars_);
	ubs_.resize(num_vars_);

	lbs_[0]	= 0.;	ubs_[0]	= 1.;	// a
	lbs_[1]	= -3.;	ubs_[1]	= 3.;	// mu1
	lbs_[2]	= -2.;	ubs_[2]	= 2.;	// mu2
	lbs_[3]	= 0.;	ubs_[3]	= 1.5;	// sigma11
	lbs_[4]	= 0.;	ubs_[4]	= 1.;	// sigma12
	lbs_[5]	= 0.;	ubs_[5]	= 1.5;	// sigma21
	lbs_[6]	= 0.;	ubs_[6]	= 1.;	// sigma22
	lbs_[7]	= -1.;	ubs_[7]	= 1.;	// kappa1
	lbs_[8]	= -1.;	ubs_[8]	= 1.;	// kappa2
	lbs_[9]	= -1.;	ubs_[9]	= 1.;	// kappa3
}
// -----------------------------------------------------------
bool CProblemSelf::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_vars_) return false;

	f.resize(num_objectives(), 0);

	/***** ここから編集 *****/
    std::size_t i;
	double beta2	= 1. / lambda_;
	// パルス振幅（generalized Gauss distribution）に関するパラメータ
	double ggd_a = sqrt(gsl_sf_gamma(1. / GGD_KAPPA)*pow(gsl_sf_gamma(3. / GGD_KAPPA), -1.)*beta2);

	// 入力に関するモーメント
	std::vector<double> dF(6);
	dF[0] = 0;
	dF[1] = alpha_*S0 + lambda_*(1. - alpha_)*beta2;
	dF[2] = 0;
	dF[3] = lambda_*pow((1. - alpha_), 2.)*(pow(ggd_a, 4.)*gsl_sf_gamma(5. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));
	dF[4] = 0;
	dF[5] = lambda_*pow((1. - alpha_), 3.)*(pow(ggd_a, 6.)*gsl_sf_gamma(7. / GGD_KAPPA)*pow(gsl_sf_gamma(1. / GGD_KAPPA), -1.));

	/******* 考慮する目的関数 *******/
	std::vector<std::size_t> rObj{0, 4, 11, 12, 13};
	// std::vector<unsigned int> rObj{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
	/******* 考慮する目的関数 *******/
    for (i = 0; i < rObj.size(); ++i) {
        f[i] = MomentEq::expb_f(x, dF, rObj)[i];
    }

	/***** ここまで編集 *****/

	return true;
}
