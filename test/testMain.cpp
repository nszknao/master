#include "../include/analysis.h"
#include <typeinfo>

int main(int argc, char *argv[])
{
	unsigned int i;
	char *ends;
	double lambda	= strtod(argv[1],&ends);
	double beta2	= strtod(argv[2],&ends);
	double alpha	= strtod(argv[3],&ends);
	double mu1		= strtod(argv[4],&ends);

	std::cout << "--------------------\n" << std::endl;
	std::cout << "analysis.cpp started.\n" << std::endl;
	
	// Parameter* prm	= new Parameter();
	// prm->allocParameter();
	Analysis *ana	= new Analysis(lambda, beta2, alpha, mu1, mu1);
	// std::string result	= ana->leastSquareMethod(prm);
	std::vector<Parameter*> prm;
	ana->GeneticAlgorithm(prm);
	// if (result == "success") {
	// 	int xmin	= -6;
	// 	std::vector<double> dispX(abs(xmin)*2*100), dispY(abs(xmin)*2*100);
	// 	ana->createDispPdf(prm, dispX, dispY, xmin);
	// 	ana->outputIntoFile((char*)"gsay1pdf.dat", dispX, dispY);
	// }
	int xmin	= -6;
	std::string filename	= "";
	for (i = 0; i < prm.size(); ++i) {
		filename	= "gsay1pdf_" + std::to_string(i) + ".dat";
		std::vector<double> dispX(abs(xmin)*2*100), dispY(abs(xmin)*2*100);
		if (!prm[i]->validate()) continue;
		ana->createDispPdf(prm[i], dispX, dispY, xmin);
		ana->outputIntoFile(filename, dispX, dispY);
	}
	
	// delete prm;
	delete ana;
	
	std::cout << "analysis.cpp has done.\n" << std::endl;
	std::cout << "--------------------\n" << std::endl;
	return 0;
}
