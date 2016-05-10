#include "../ga/GA.h"


int main(int argc, char *argv[])
{

        int generation = 0, numVarience = 10;

        GA ga(numVarience);
        ga->nsga2Run();

        return 0;
}