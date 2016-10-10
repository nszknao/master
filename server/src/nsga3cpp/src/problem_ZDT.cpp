#include "problem_ZDT.h"
#include "alg_individual.h"
#include "aux_math.h"
#include <cmath>
#include <numeric>

using namespace std;


// ----------------------------------------------------------------------
//		CProblemZDT1
// ----------------------------------------------------------------------

CProblemZDT1::CProblemZDT1(std::size_t n):
	CProblemZDT(n, "ZDT1")
{
	lbs_.resize(num_variables_, 0.0); // lower bound 0.0
	ubs_.resize(num_variables_, 1.0); // upper bound 1.0
}
// ----------------------------------------------------------------------
bool CProblemZDT1::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_variables_) return false; // #variables does not match

	f.resize(num_objectives(), 0);

	f[0] = x[0];

	double g = accumulate(x.begin()+1, x.end(), 0.0);
	g = 1 + 9*g/(num_variables_ - 1);

	double h = 1 - sqrt(x[0]/g);
	f[1] = g*h;

	return true;

}// bool CProblemZDT1::Evaluate(CIndividual *indv) const



// ----------------------------------------------------------------------
//		CProblemZDT2
// ----------------------------------------------------------------------

CProblemZDT2::CProblemZDT2(std::size_t n):
	CProblemZDT(n, "ZDT2")
{
	lbs_.resize(num_variables_, 0.0); // lower bound 0.0
	ubs_.resize(num_variables_, 1.0); // upper bound 1.0
}
// ----------------------------------------------------------------------
bool CProblemZDT2::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_variables_) return false; // #variables does not match

	f.resize(num_objectives(), 0);

	f[0] = x[0];

	double g = accumulate(x.begin()+1, x.end(), 0.0);
	g = 1 + 9*g/(num_variables_ - 1);

	double h = 1 - MathAux::square(x[0]/g);
	f[1] = g*h;

	return true;

}// bool CProblemZDT2::Evaluate(CIndividual *indv) const


// ----------------------------------------------------------------------
//		CProblemZDT3
// ----------------------------------------------------------------------

CProblemZDT3::CProblemZDT3(std::size_t n):
	CProblemZDT(n, "ZDT3")
{
	lbs_.resize(num_variables_, 0.0); // lower bound 0.0
	ubs_.resize(num_variables_, 1.0); // upper bound 1.0
}
// ----------------------------------------------------------------------
bool CProblemZDT3::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_variables_) return false; // #variables does not match

	f.resize(num_objectives(), 0);

	f[0] = x[0];

	double g = accumulate(x.begin()+1, x.end(), 0.0);
	g = 1 + 9*g/(num_variables_ - 1);

	double h = 1 - sqrt(x[0]/g) - x[0]/g*sin(10*MathAux::PI*x[0]);
	f[1] = g*h;

	return true;

}// bool CProblemZDT3::Evaluate(CIndividual *indv) const



// ----------------------------------------------------------------------
//		CProblemZDT4
// ----------------------------------------------------------------------

CProblemZDT4::CProblemZDT4(std::size_t n):
	CProblemZDT(n, "ZDT4")
{
	lbs_.resize(num_variables_, -5.0); // lower bound
	ubs_.resize(num_variables_, 5.0); // upper bound

	lbs_[0] = 0.0;
	ubs_[0] = 1.0;
}
// ----------------------------------------------------------------------
bool CProblemZDT4::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_variables_) return false; // #variables does not match

	f.resize(num_objectives(), 0);

	f[0] = x[0];

	double g = 1 + 10*(num_variables_ - 1);
	
	for (size_t i=1; i<num_variables_; i+=1)
	{
		g += x[i]*x[i] - 10*cos(4*MathAux::PI*x[i]);
	}

	double h = 1 - sqrt(x[0]/g);
	f[1] = g*h;

	return true;

}// bool CProblemZDT4::Evaluate(CIndividual *indv) const


// ----------------------------------------------------------------------
//		CProblemZDT6
// ----------------------------------------------------------------------

CProblemZDT6::CProblemZDT6(std::size_t n):
	CProblemZDT(n, "ZDT6")
{
	lbs_.resize(num_variables_, 0.0); // lower bound
	ubs_.resize(num_variables_, 1.0); // upper bound
}
// ----------------------------------------------------------------------
bool CProblemZDT6::Evaluate(CIndividual *indv) const
{
	CIndividual::TDecVec &x = indv->vars();
	CIndividual::TObjVec &f = indv->objs();

	if (x.size() != num_variables_) return false; // #variables does not match

	f.resize(num_objectives(), 0);

	f[0] = 1 - exp(-4*x[0])*pow(sin(4*MathAux::PI*x[0]), 6);

	double g = accumulate(x.begin()+1, x.end(), 0.0);
	g = 1 + 9*pow(g/(num_variables_ - 1), 0.25);

	double h = 1 - MathAux::square(f[0]/g);
	f[1] = g*h;

	return true;

}// bool CProblemZDT6::Evaluate(CIndividual *indv) const