#ifndef PROBLEM_SELF__
#define PROBLEM_SELF__

#include "problem_base.h"
#include <gsl/gsl_sf.h>

class MomentEq;

class CProblemSelf : public BProblem
{
public:
	CProblemSelf(MomentEq *);

	virtual std::size_t num_variables() const { return num_vars_; }
	virtual std::size_t num_objectives() const { return num_objs_; }
	virtual bool Evaluate(CIndividual *indv) const;

private:
    std::size_t num_vars_, num_objs_;
    MomentEq *meq_;
};

#endif
