
#include "problem_factory.h"
#include "problem_DTLZ.h"
#include "problem_ZDT.h"
#include "problem_self.h"
#include <string>
using namespace std;

BProblem *GenerateProblem(MomentEq *meq)
{
	return new CProblemSelf(meq);
}
