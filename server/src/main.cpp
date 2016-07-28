#include "../include/analysis.h"
#include "../include/research.h"
#include "../include/common.h"


int main(int argc, char *argv[])
{
	unsigned int i;
	std::string filename	= "";

	char *ends;
	double lambda	= strtod(argv[1],&ends);
	double beta2	= strtod(argv[2],&ends);
	double alpha	= strtod(argv[3],&ends);
	double mu1		= strtod(argv[4],&ends);

	std::cout << "--------------------\n" << std::endl;
	std::cout << "research.cpp started.\n" << std::endl;

	Simulation * sim	= new Simulation(lambda, beta2, alpha);

	/* ルンゲクッタを解く */
	std::vector<double> force, t;
	std::vector< std::vector<double> > y1, y2;
	sim->culcRungeKutta(t, y1, y2, force);
	filename	= "t_force.dat";
	Common::outputIntoFile(filename, t, force);
	filename	= "t_x1.dat";
	Common::outputIntoFile(filename, t, y1[0]);
	/* 変位のPDFを求める */
	std::vector<double> dispX, dispY;
	sim->createDispPdf(y1, dispX, dispY);
	filename	= "y1_pdf.dat";
	Common::outputIntoFile(filename, dispX, dispY);
	delete sim;

	std::cout << "analysis.cpp has done.\n" << std::endl;
	std::cout << "--------------------\n" << std::endl;


	std::cout << "--------------------\n" << std::endl;
	std::cout << "analysis.cpp started.\n" << std::endl;
	
	Analysis *ana	= new Analysis(lambda, beta2, alpha, mu1, mu1);

	/* 最小二乗法で解く */
	Parameter* prm	= new Parameter();
	prm->allocParameter();
	std::string result	= ana->leastSquareMethod(prm);
	if (result == "success") {
		int xmin	= -6;
		std::vector<double> dispX(abs(xmin)*2*100), dispY(abs(xmin)*2*100);
		ana->createDispPdf(prm, dispX, dispY, xmin);
		ana->outputIntoFile((char*)"gsay1pdf.dat", dispX, dispY);
	}
	delete prm;

	/* GAで解く */
	// std::vector<Parameter*> prm;
	// std::vector< std::vector<double> > pValue, oValue;
	// ana->GeneticAlgorithm(prm, pValue, oValue);
	// int xminDisp	= -6;
	// int xmaxFCross	= 8;
	// for (i = 0; i < prm.size(); ++i) {
	// 	if (!prm[i]->validate()) continue;
	// 	// if (ana->isOverSpecifyValue(oValue[i], 50.)) continue;
	// 	/* 変位のPDFを求める */
	// 	std::vector<double> dispX(abs(xminDisp)*2*100), dispY(abs(xminDisp)*2*100);
	// 	ana->createDispPdf(prm[i], dispX, dispY, xminDisp);
	// 	filename	= "gsay1pdf_" + std::to_string(i) + ".dat";
	// 	Common::outputIntoFile(filename, dispX, dispY);
	// 	/* 閾値通過率を求める */
	// 	std::vector<double> fCrossX(abs(xmaxFCross)*100), fCrossY(abs(xmaxFCross)*100);
	// 	ana->createLevelCrossing(prm[i], fCrossX, fCrossY, xmaxFCross);
	// 	filename	= "firstcross_" + std::to_string(i) + ".dat";
	// 	// Common::outputIntoFile(filename, fCrossX, fCrossY);
	// }
	// for (i = 0; i < prm.size(); ++i) {
	// 	prm[i]->freeParameter();
	// 	delete prm[i];
	// }

	delete ana;
	
	std::cout << "analysis.cpp has done.\n" << std::endl;
	std::cout << "--------------------\n" << std::endl;
	return 0;
}
