#ifndef __ProblemSelf_H__
#define __ProblemSelf_H__

#include <Problem.h>
#include <math.h>
#include <BinaryRealSolutionType.h>
#include <RealSolutionType.h>
#include <ArrayRealSolutionType.h>
#include <XReal.h>
#include <Solution.h>

#include <expfit.h>
#include <common.h>

class ProblemSelf : public Problem {
public:
    ProblemSelf(string solutionType, MomentEq &);
    void evaluate(Solution *solution);

    virtual ~ProblemSelf();
private:
    MomentEq meq_;
//    double * fx_ ;
//    double * x_  ;
    std::vector<double> fx_;
    std::vector<double> x_;
};

#endif /* __ProblemSelf_H__ */
