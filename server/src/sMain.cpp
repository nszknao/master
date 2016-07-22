#include "../include/research.h"
#include "../include/common.h"


int main(int argc, char *argv[])
{
	std::string filename	= "";

	char *ends;
	double lambda	= strtod(argv[1],&ends);
	double beta2	= strtod(argv[2],&ends);
	double alpha	= strtod(argv[3],&ends);

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

	return 0;
}
