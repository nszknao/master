
#ifndef EXPERIMENT__
#define EXPERIMENT__

#include <fstream>
#include "../../../include/expfit.h"

class CNSGAIII;
class BProblem;

void SetupExperiment(CNSGAIII &algo, BProblem **prob, std::ifstream &ifile, MomentEq *meq);

#endif
