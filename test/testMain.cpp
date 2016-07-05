#include "../include/analysis.h"


int main(int argc, char *argv[])
{
		char *ends;
		double lambda	= strtod(argv[1],&ends);
		double beta2	= strtod(argv[2],&ends);
		double alpha	= strtod(argv[3],&ends);
		double mu1		= strtod(argv[4],&ends);
        Analysis *ana	= new Analysis(lambda, beta2, alpha, mu1, mu1);

        std::cout << "--------------------\n" << std::endl;
        std::cout << "analysis.cpp started.\n" << std::endl;

        Parameter* prm	= new Parameter();
        prm->allocParameter();
        std::string result	= ana->leastSquareMethod(prm);

        if (result == "success") {
	       	ana->createDispPdf(prm);
        }

        prm->freeParameter();
        delete prm;

        std::cout << "analysis.cpp has done.\n" << std::endl;
        std::cout << "--------------------\n" << std::endl;

        return 0;
}
