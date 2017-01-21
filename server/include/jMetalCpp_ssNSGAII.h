#ifndef __jMetalCpp_ssNSGAII_H_INCLUDE_
#define __jMetalCpp_ssNSGAII_H_INCLUDE_

#include <Problem.h>
#include <Solution.h>
#include <SBXCrossover.h>
#include <PolynomialMutation.h>
#include <BinaryTournament2.h>
#include <iostream>
#include <NSGAII.h>
#include <ssNSGAII.h>
#include <ProblemFactory.h>
#include <string.h>
#include <time.h>

#include <ga_individual.h>
#include <expfit.h>

using namespace std;

class GAIndividual;
class MomentEq;

class jMetalCpp_ssNSGAII
{
private:
    std::vector<GAIndividual> _finalPops;
    std::size_t _popSize, _iterations;

    void _saveArchive(SolutionSet &);

public:
    explicit jMetalCpp_ssNSGAII(std::size_t, std::size_t);
    int run(MomentEq *);
    std::vector<GAIndividual> getFinalPops();
    ~jMetalCpp_ssNSGAII();
};

#endif // !__jMetalCpp_ssNSGAII_H_INCLUDE_
