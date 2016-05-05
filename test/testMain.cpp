#include "../ga/GA.h"


int main(int argc, char *argv[])
{

        int generation = 0, numVarience = 1;
        double mutationRate = 0.03;

        GA ga(numVarience);

        ga.initGene();
        ga.culcFitness();

        ga.outputGeneration(generation);

        for (generation = 1; generation < 21; generation++)
        {
                ga.selectRanking();
                ga.uniformCrossover();
                ga.mutation(mutationRate);
                ga.culcFitness();

                ga.outputGeneration(generation);
        }

        return 0;
}