#ifndef __NSGA3_H_INCLUDE_
#define __NSGA3_H_INCLUDE_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class CPopulation;
class GAIndividual;
class MomentEq;

class NSGA3
{
private:
    std::vector<GAIndividual> _finalPops;
    void _savePopulation(const CPopulation &);

public:
    ~NSGA3();
    int run(MomentEq *);
    std::vector<GAIndividual> getFinalPops();
};

#endif // !__NSGA3_H_INCLUDE_
