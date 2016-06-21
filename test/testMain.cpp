#include "../include/analysis.h"


int main(int argc, char *argv[])
{
		char *ends;
        Analysis *ana;
        ana	= new Analysis(strtod(argv[1], &ends),strtod(argv[2], &ends),strtod(argv[3], &ends),strtod(argv[4], &ends),strtod(argv[4], &ends));

        std::cout << "--------------------\n" << std::endl;
        std::cout << "analysis.cpp started.\n" << std::endl;

        std::string result	= ana->leastSquareMethod();
        if (result == "success") {
        	ana->createDispPdf();
        }

        std::cout << "analysis.cpp has done.\n" << std::endl;
        std::cout << "--------------------\n" << std::endl;

        return 0;
}
