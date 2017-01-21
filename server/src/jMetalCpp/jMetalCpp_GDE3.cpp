//  GDE3_main.cpp
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

#include <jMetalCpp_GDE3.h>

jMetalCpp_GDE3::jMetalCpp_GDE3(std::size_t pop, std::size_t iter)
:_popSize(pop),
_iterations(iter)
{}

jMetalCpp_GDE3::~jMetalCpp_GDE3()
{
      std::vector<GAIndividual>().swap(_finalPops);
}

int jMetalCpp_GDE3::run(MomentEq *meq) {

    clock_t t_ini, t_fin;

    Problem   * problem   ; // The problem to solve
    Algorithm * algorithm ; // The algorithm to use
    Operator  * crossover ; // Crossover operator
    Operator  * selection ; // Selection operator

    map<string, void *> parameters;
    
    //TODO: QualityIndicator * indicators;

    // 引数指定しなければexpfitが呼ばれるようにした
    problem = ProblemFactory::getProblem(*meq);
    cout << "Selected problem: " << problem->getName() << endl;
 
    algorithm = new GDE3(problem);

    // Algorithm parameters
    int populationSizeValue = _popSize;
    int maxIterationsValue = _iterations;
    algorithm->setInputParameter("populationSize",&populationSizeValue);
    algorithm->setInputParameter("maxIterations",&maxIterationsValue);

    // Crossover operator
    double crParameter = 0.5;
    double fParameter  = 0.5;
    parameters["CR"] =  &crParameter;
    parameters["F"] = &fParameter;
    crossover = new DifferentialEvolutionCrossover(parameters);

    // Selection operator
    parameters.clear();
    selection = new DifferentialEvolutionSelection(parameters) ;

    // Add the operators to the algorithm
    algorithm->addOperator("crossover",crossover);
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
    cout << "Variables values have been written to file VAR" << endl;
    population->printVariablesToFile("VAR");
    cout << "Objectives values have been written to file FUN" << endl;
    population->printObjectivesToFile("FUN");

    delete selection;
    delete crossover;
    delete population;
    delete algorithm;

} // main

/**
 * @fn 解の個体群を取得
 * @return vector<GAIndividual> _finalPops 最終的な個体群
 */
std::vector<GAIndividual> jMetalCpp_GDE3::getFinalPops()
{
    return _finalPops;
}

/**
 * @fn アーカイブの情報を_finalPopsへ格納
 * @param ArchiveMOO &archive アーカイブ集団
 */
void jMetalCpp_GDE3::_saveArchive(SolutionSet &pop)
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
