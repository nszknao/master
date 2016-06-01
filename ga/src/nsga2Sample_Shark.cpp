#include <stdlib.h>
#include <stdio.h>

#include <Array/Array.h>
#include <FileUtil/Params.h>
#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>

using namespace std;


inline void calcFitness(IndividualMOO &i)
{
    double f1, f2;
//  ChromosomeT< double >& dblchrom = dynamic_cast< ChromosomeT< double >& >( i[ 0 ] );

    f1 = Rng::uni(0, 1);
    f2 = Rng::uni(0, 1);

    i.setMOOFitnessValues(f1, f2);
}

// main program
int main(int argc, char* argv[])
{
    unsigned i, ii, t;

    unsigned seed        = 0;
    unsigned dimension   = 10;  // problem dimension
    unsigned mu          = 100; // population size
    unsigned iterations  = 100; // maximum number of generations
    unsigned archiveSize = mu;  // how many solutions in external archive

    // bounds
    double   minInit = -1;
    double   maxInit = 1;

    // parameters of NSGA2 mutation and crossover operators
    double   nm = 20;
    double   nc = 20;
    double   pc = 0.9;            // crossover probability
    double   pm = 1. / dimension; // mutation probability

    // initialize random number generator
    Rng::seed(seed);

    // define populations
    PopulationMOO parents(mu, ChromosomeT< double >(dimension));
    PopulationMOO offsprings(mu, ChromosomeT< double >(dimension));

    // set minimization task
    parents   .setMinimize();
    offsprings.setMinimize();

    // set # objective funstions
    parents   .setNoOfObj(2);
    offsprings.setNoOfObj(2);

    // init archive
    ArchiveMOO archive(archiveSize);
    archive.minimize();
    archive.setMaxArchive(archiveSize);

    // bounds used by operators
    vector<double> lower(dimension), upper(dimension);
    for (i = 0; i < lower.size(); ++i) {
        lower[i] = minInit;
        upper[i] = maxInit;
    }

    // initialize all chromosomes of parent population
    for (i = 0; i < parents.size(); ++i)
//    dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize( minInit, maxInit);
        ((ChromosomeT< double >*)&parents[ i ][ 0 ])->initialize(minInit, maxInit);

    // evaluate parents
    for (i = 0; i < parents.size(); ++i)
        calcFitness(parents[ i ]);

    parents.crowdedDistance();

    // iterate
    for (t = 1; t <= iterations; ++t) {
        cout << "generation: " << t << endl;

        // copy parents to offsprings
        offsprings.selectBinaryTournamentMOO(parents);

        // recombine by crossing over two parents
        for (i = 0; i < offsprings.size(); i += 2)
            if (Rng::coinToss(pc))
                dynamic_cast< ChromosomeT< double >& >(offsprings[i][0]).SBX(dynamic_cast< ChromosomeT< double >& >(offsprings[i+1][0]), lower, upper, nc, .5);

        // mutate by flipping bits and evaluate objective function
        for (i = 0; i < offsprings.size(); ++i) {
            // modify
            dynamic_cast< ChromosomeT< double >& >(offsprings[i][0]).mutatePolynomial(lower, upper, nm, pm);

            // evaluate objective function
            calcFitness(offsprings[ i ]);
        }

        // selection
        parents.selectCrowdedMuPlusLambda(offsprings);

        // every 10 generations the dominated solutions are removed from the archive
        if (!(t % 10)) {
            archive.cleanArchive();
            for (i = 0; i < parents.size(); i++)
                archive.addArchive(parents[i]);
            archive.nonDominatedSolutions();
        }
    }

    // data output
    archive.cleanArchive();
    for (ii = 0; ii < archiveSize; ii++)
        archive.addArchive(parents[ ii ]);

    archive.nonDominatedSolutions();

    char filename[256];
    sprintf(filename, "nsga2example-%d.out", seed);
    cout    << "size of the archive: "  << archive.size()
    << ", filename of archive: " << filename << endl << endl;

    // ファイルに目的関数値を書き込む
    archive.saveArchiveGPT(filename);

    return EXIT_SUCCESS;
}
