
#include "exp_experiment.h"
#include "problem_factory.h"
#include "alg_nsgaiii.h"
#include "../../../include/expfit.h"

void SetupExperiment(CNSGAIII &algo, BProblem **prob, std::ifstream &ifile, MomentEq *meq)
{
	algo.Setup(ifile);
	*prob = GenerateProblem(meq);
}
