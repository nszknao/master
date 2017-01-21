#ifndef __jMetalCpp_GDE3_H_INCLUDE_
#define __jMetalCpp_GDE3_H_INCLUDE_

#include <Problem.h>
#include <Solution.h>
#include <SBXCrossover.h>
#include <PolynomialMutation.h>
#include <BinaryTournament2.h>
#include <DifferentialEvolutionCrossover.h>
#include <DifferentialEvolutionSelection.h>
#include <iostream>
#include <GDE3.h>
#include <ProblemFactory.h>
#include <string.h>
#include <time.h>

#include <ga_individual.h>
#include <expfit.h>

using namespace std;

class GAIndividual;
class MomentEq;

class jMetalCpp_GDE3
{
private:
    std::vector<GAIndividual> _finalPops;
    std::size_t _popSize, _iterations;

    void _saveArchive(SolutionSet &);

public:
    explicit jMetalCpp_GDE3(std::size_t, std::size_t);
    int run(MomentEq *);
    std::vector<GAIndividual> getFinalPops();
    ~jMetalCpp_GDE3();
};

#endif // !__jMetalCpp_GDE3_H_INCLUDE_
