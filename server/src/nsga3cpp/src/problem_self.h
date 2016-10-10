#ifndef PROBLEM_SELF__
#define PROBLEM_SELF__

#include "problem_base.h"
#include <gsl/gsl_sf.h>

/********** 系の係数・入力条件**********/
#ifndef __OSCILLATOR_AND_EXCITATION_PARAMETER__
#define __OSCILLATOR_AND_EXCITATION_PARAMETER__
#define PI M_PI
#define S0 1./(2.*PI)
#define EPSILON 0.3
#define ZETA 0.05
#endif // !__OSCILLATOR_AND_EXCITATION_PARAMETER__

#define GGD_KAPPA 2.	// 1.:ラプラス分布，2.:ガウス分布，∞:一様分布

class CProblemSelf : public BProblem
{
public:
	CProblemSelf(std::size_t num_vars, std::size_t num_objs, double lambda, double alpha);

	virtual std::size_t num_variables() const { return num_vars_; }
	virtual std::size_t num_objectives() const { return num_objs_; }
	virtual bool Evaluate(CIndividual *indv) const;

private:

	std::size_t num_vars_, num_objs_;
	double lambda_, alpha_;
};

#endif