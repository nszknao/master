#ifndef PROBLEM_ZDT__
#define PROBLEM_ZDT__

#include "problem_base.h"
#include <cstddef>


// ----------------------------------------------------------------------
//		CProblemZDT (2-objective problems)
//
// E. Zitzler, K. Deb, and L. Thiele, "Comparison of Multiobjective 
// Evolutionary Algorithms: Empirical Results," Evolutionary Computation,
// vol. 8, no. 2, pp. 173-195, 2000.
// ----------------------------------------------------------------------
class CProblemZDT : public BProblem
{
public:
	explicit CProblemZDT(std::size_t n, const std::string &name):
	BProblem(name),
	num_variables_(n)
	{}

	virtual std::size_t num_variables() const { return num_variables_; }
	virtual std::size_t num_objectives() const { return 2; }

protected:
	std::size_t num_variables_;
};


// ----------------------------------------------------------------------
//		CProblemZDT1
// ----------------------------------------------------------------------
class CProblemZDT1 : public CProblemZDT
{
public:
	explicit CProblemZDT1(std::size_t n = 30); // default setting in ZDT1
	virtual bool Evaluate(CIndividual *indv) const;
};

// ----------------------------------------------------------------------
//		CProblemZDT2
// ----------------------------------------------------------------------
class CProblemZDT2 : public CProblemZDT
{
public:
	explicit CProblemZDT2(std::size_t n = 30); // default setting in ZDT2
	virtual bool Evaluate(CIndividual *indv) const;
};

// ----------------------------------------------------------------------
//		CProblemZDT3
// ----------------------------------------------------------------------
class CProblemZDT3 : public CProblemZDT
{
public:
	explicit CProblemZDT3(std::size_t n = 30); // default setting in ZDT3
	virtual bool Evaluate(CIndividual *indv) const;
};

// ----------------------------------------------------------------------
//		CProblemZDT4
// ----------------------------------------------------------------------
class CProblemZDT4 : public CProblemZDT
{
public:
	explicit CProblemZDT4(std::size_t n = 10); // default setting in ZDT4
	virtual bool Evaluate(CIndividual *indv) const;
};


// ----------------------------------------------------------------------
//		CProblemZDT6
// ----------------------------------------------------------------------
class CProblemZDT6 : public CProblemZDT
{
public:
	explicit CProblemZDT6(std::size_t n = 10); // default setting in ZDT6
	virtual bool Evaluate(CIndividual *indv) const;
};

#endif