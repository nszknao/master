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
    lbs_.resize(Common::NUM_OF_PARAM);
    ubs_.resize(Common::NUM_OF_PARAM);

    lbs_[0] = 0.3;   ubs_[0] = 1.;   // a
    lbs_[1] = -3.;  ubs_[1] = 3.;   // mu1
    lbs_[2] = -2.;  ubs_[2] = 2.;   // mu2
    lbs_[3] = 0.;   ubs_[3] = 1.5;  // sigma11
    lbs_[4] = 0.;   ubs_[4] = 1.;   // sigma12
    lbs_[5] = 0.;   ubs_[5] = 1.5;  // sigma21
    lbs_[6] = 0.;   ubs_[6] = 1.;   // sigma22
    lbs_[7] = -1.;  ubs_[7] = 1.;   // rho1
    lbs_[8] = -1.;  ubs_[8] = 1.;   // rho2
    lbs_[9] = -1.;  ubs_[9] = 1.;   // rho3
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
