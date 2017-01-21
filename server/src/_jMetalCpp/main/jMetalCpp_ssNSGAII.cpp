//  ssNSGAII_main.cpp
//
//  Author:
//       Esteban López-Camacho <esteban@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <jMetalCpp_ssNSGAII.h>

jMetalCpp_ssNSGAII::jMetalCpp_ssNSGAII(std::size_t pop, std::size_t iter)
:_popSize(pop),
_iterations(iter)
{}

jMetalCpp_ssNSGAII::~jMetalCpp_ssNSGAII()
{
    std::vector<GAIndividual>().swap(_finalPops);
}

/**
 * Class implementing the steady-state version of the NSGA-II algorithm.
 * This implementation of ssNSGA-II makes use of a QualityIndicator object
 *  to obtained the convergence speed of the algorithm. This version is used
 *  in the paper:
 *     J.J. Durillo, A.J. Nebro, F. Luna and E. Alba
 *     "On the Effect of the Steady-State Selection Scheme in Multi-Objective
 *     Genetic Algorithms" 5th International Conference, EMO 2009, pp: 183-197.
 *     April 2009)
 */
int jMetalCpp_ssNSGAII::run(MomentEq *meq) {

    clock_t t_ini, t_fin;
    
    Problem   * problem   ; // The problem to solve
    Algorithm * algorithm ; // The algorithm to use
    Operator  * crossover ; // Crossover operator
    Operator  * mutation  ; // Mutation operator
    Operator  * selection ; // Selection operator

    // 引数指定しなければexpfitが呼ばれるようにした
    problem = ProblemFactory::getProblem(*meq);
    cout << "Selected problem: " << problem->getName() << endl;
 
   //TODO: Quality Indicators
    //QualityIndicator indicators ; // Object to get quality indicators
    //indicators = null ;

    algorithm = new ssNSGAII(problem);

    // Algorithm parameters
    int populationSize = _popSize;
    int maxEvaluations = _iterations;
    algorithm->setInputParameter("populationSize",&populationSize);
    algorithm->setInputParameter("maxEvaluations",&maxEvaluations);

    // Mutation and Crossover for Real codification
    map<string, void *> parameters;

    double crossoverProbability = 0.9;
    double crossoverDistributionIndex = 20.0;
    parameters["probability"] =  &crossoverProbability;
    parameters["distributionIndex"] = &crossoverDistributionIndex;
    crossover = new SBXCrossover(parameters);

    parameters.clear();
//    double mutationProbability = 1.0/problem->getNumberOfVariables();
    double mutationProbability = 0.01;
    double mutationDistributionIndex = 20.0;
    parameters["probability"] = &mutationProbability;
    parameters["distributionIndex"] = &mutationDistributionIndex;
    mutation = new PolynomialMutation(parameters);

    // Selection Operator
    parameters.clear();
    selection = new BinaryTournament2(parameters);

    // Add the operators to the algorithm
    algorithm->addOperator("crossover",crossover);
    algorithm->addOperator("mutation",mutation);
    algorithm->addOperator("selection",selection);

    // Add the indicator object to the algorithm
    //algorithm->setInputParameter("indicators", indicators) ;

    // Execute the Algorithm
    t_ini = clock();
    SolutionSet * population = algorithm->execute();
    t_fin = clock();
    double secs = (double) (t_fin - t_ini);
    secs = secs / CLOCKS_PER_SEC;
    _saveArchive(*population);

    // Result messages
    cout << "Total execution time: " << secs << "s" << endl;
//    cout << "Variables values have been written to file VAR" << endl;
//    population->printVariablesToFile("VAR");
//    cout << "Objectives values have been written to file FUN" << endl;
//    population->printObjectivesToFile("FUN");

    delete selection;
    delete mutation;
    delete crossover;
    delete population;
    delete algorithm;

    return EXIT_SUCCESS;
} // main

/**
 * @fn 解の個体群を取得
 * @return vector<GAIndividual> _finalPops 最終的な個体群
 */
std::vector<GAIndividual> jMetalCpp_ssNSGAII::getFinalPops()
{
    return _finalPops;
}

/**
 * @fn アーカイブの情報を_finalPopsへ格納
 * @param ArchiveMOO &archive アーカイブ集団
 */
void jMetalCpp_ssNSGAII::_saveArchive(SolutionSet &pop)
{
    std::size_t i, ii;
    std::size_t no = pop.size();
    std::size_t noOfObj;
    if (no > 0)
        noOfObj = pop.get(0)->getNumberOfObjectives();
    else
        noOfObj = 0;

    _finalPops.resize(no);
    for (i = 0; i < no; ++i) {
        std::vector<double> variable(pop.get(i)->getNumberOfVariables());
        for (ii = 0; ii < variable.size(); ++ii) {
            variable[ii] = pop.get(i)->getDecisionVariables()[ii]->getValue();
        }
        // モーメント値
        std::vector<double> m = MomentEq::getMomentFromParameter(variable);
        _finalPops[i].mValue.resize(m.size());
        for (ii = 0; ii < m.size(); ++ii) {
            _finalPops[i].mValue[ii]    = m[ii];
        }
        // 目的関数値
        _finalPops[i].oValue.resize(noOfObj);
        for (ii = 0; ii < noOfObj; ++ii) {
            _finalPops[i].oValue[ii] = pop.get(i)->getObjective(ii);
        }
        // パラメータ値
        _finalPops[i].pValue.resize(variable.size());
        for (ii = 0; ii < variable.size(); ++ii) {
            _finalPops[i].pValue[ii] = variable[ii];
        }
    }
}
