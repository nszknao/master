#include "problem_self.h"
#include "alg_individual.h"

#include "../../../include/expfit.h"
#include "../../../include/common.h"

using namespace std;

CProblemSelf::CProblemSelf(MomentEq *meq):
    BProblem("Self"),
    num_vars_(Common::NUM_OF_PARAM),
    num_objs_(meq->getObjNum()),
    meq_(meq)
{
    // define the domain of variables here
    lbs_ = meq_->getLowerObj();
    ubs_ = meq_->getUpperObj();
}
// -----------------------------------------------------------
bool CProblemSelf::Evaluate(CIndividual *indv) const
{
    CIndividual::TDecVec &x = indv->vars();
    CIndividual::TObjVec &f = indv->objs();

    if (x.size() != Common::NUM_OF_PARAM) return false;

    f.resize(meq_->getObjNum(), 0);

    std::size_t i;
    for (i = 0; i < meq_->getObjNum(); ++i) {
        f[i] = meq_->expb_f(x)[i];
    }

    return true;
}
