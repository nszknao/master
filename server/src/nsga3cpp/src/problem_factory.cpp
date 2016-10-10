
#include "problem_factory.h"
#include "problem_DTLZ.h"
#include "problem_ZDT.h"
#include "problem_self.h"
#include <string>
using namespace std;

BProblem *GenerateProblem(std::ifstream &ifile)
{

	string pname, dummy;

	ifile >> dummy >> dummy >> pname; // proble.name = [ ]
	if (pname == "Self")
	{
		// define your own file format, read it, and generate an instance here

		size_t num_vars, num_objs;
		double lambda, alpha;

		ifile >> dummy >> dummy >> num_vars;
		ifile >> dummy >> dummy >> num_objs;
		ifile >> dummy >> dummy >> lambda;
		ifile >> dummy >> dummy >> alpha;

		return new CProblemSelf(num_vars, num_objs, lambda, alpha);
	}
	else if (pname.substr(0, 4) == "DTLZ")
	{
		size_t num_objs;

		ifile >> dummy >> dummy >> num_objs; // problem.objectives = [ ]

		if (pname == "DTLZ1") return new CProblemDTLZ1(num_objs);
		else if (pname == "DTLZ2") return new CProblemDTLZ2(num_objs);
		else if (pname == "DTLZ3") return new CProblemDTLZ3(num_objs);
		else if (pname == "DTLZ4") return new CProblemDTLZ4(num_objs);
		else return 0;
	}
	else if (pname.substr(0, 3) == "ZDT")
	{
		size_t num_objs;

		ifile >> dummy >> dummy >> num_objs; // problem.objectives = [ ]

		if (pname == "ZDT1") return new CProblemZDT1();
		else if (pname == "ZDT2") return new CProblemZDT2();
		else if (pname == "ZDT3") return new CProblemZDT3();
		else if (pname == "ZDT4") return new CProblemZDT4();
		else if (pname == "ZDT6") return new CProblemZDT6();
		else return 0;
	}
		
	return 0;
}
