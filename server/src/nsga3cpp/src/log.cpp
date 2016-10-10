
#include "log.h"
#include "alg_population.h"
#include "gnuplot_interface.h"

#include <fstream>
#include <sstream>
using namespace std;



bool SaveScatterData(const std::string &fname, const CPopulation &pop, LogDetails details, ios_base::openmode mode)
{
	ofstream ofile(fname.c_str(), mode);
	if (!ofile) return false;

	for (size_t i=0; i<pop.size(); i+=1)
	{

		if (details == ldAll)
		{
			for (size_t j=0; j<pop[i].vars().size(); j+=1)
			{
				ofile << pop[i].vars()[j] << ' ';
			}
		}

		for (size_t f=0; f<pop[i].objs().size(); f+=1)
		{
			ofile << pop[i].objs()[f] << ' ';
		}
		ofile << endl;
	}
	ofile << endl;
	return true;
}
// ----------------------------------------------------------------------
bool SaveValuePathData(const std::string &fname, const CPopulation &pop, ios_base::openmode mode)
{
	ofstream ofile(fname.c_str(), mode);
	if (!ofile) return false;

	for (size_t f=0; f<pop[0].objs().size(); f+=1)
	{
		ofile << f+1 << ' ';
		for (size_t i=0; i<pop.size(); i+=1)
		{
			ofile << pop[i].objs()[f] << ' ';
		}
		ofile << endl;
	}
	ofile << endl;

	return true;
}
// ----------------------------------------------------------------------
bool ShowPopulation(Gnuplot &gplot, const CPopulation &pop, const std::string &legend)
{
	bool successful = true;
	if (pop[0].objs().size() == 2) // 2-D scatter plot
	{
		LogDetails ld = ldObjective;
		if ((successful = SaveScatterData(legend, pop, ld)))
		{
			size_t col = (ld==ldObjective?0:pop[0].vars().size());
			gplot.plot(legend, col+1, col+2);
		}
	}
	else if (pop[0].objs().size() == 3) // 3-D scatter plot
	{
		LogDetails ld = ldObjective;
		if ((successful = SaveScatterData(legend, pop, ldObjective)))
		{
			size_t col = (ld==ldObjective?0:pop[0].vars().size());
			gplot.splot(legend, col+1, col+2, col+3);
		}
	}
	else // value path plot
	{
		gplot("set xtics 1");
		ostringstream command;
		if ((successful = SaveValuePathData(legend, pop, ios_base::out)))
		{
			command << "plot \"" << legend << "\" using 1:2 with lines title \"\", ";               
			for (size_t i=1; i<pop.size()-1; i+=1)
			{
				command << "\"\" using 1:" << i+2 << " with lines title \"\", ";
			}
			command << "\"\" using 1:" << pop.size()+1 << " with lines title \"\" ";
			gplot(command.str());
		}
		gplot("set xtics auto");
	}
	return successful;
}