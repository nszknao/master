#include <ProblemSelf.h>

ProblemSelf::ProblemSelf(string solutionType, MomentEq &meq)
:meq_(meq)
{
    std::size_t i;

    numberOfVariables_   = Common::NUM_OF_PARAM;
    numberOfObjectives_  = meq.getObjNum();
    problemName_                 = "Self";

    lowerLimit_ = new double[numberOfVariables_];
    if (lowerLimit_ == NULL) {
        cout << "Impossible to reserve memory for storing the variable lower limits" << endl;
        exit(-1);
    }

    upperLimit_ = new double[numberOfVariables_];
    if (upperLimit_ == NULL) {
        cout << "Impossible to reserve memory for storing the variable lower limits" << endl;
        exit(-1);
    }

    for (i = 0; i < numberOfVariables_; ++i) {
        lowerLimit_[i] = meq.getLowerObj()[i];
        upperLimit_[i] = meq.getUpperObj()[i];
    }

    if (solutionType.compare("BinaryReal") == 0)
        solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compare("Real") == 0) {
        solutionType_ = new RealSolutionType(this) ;
        //cout << "Tipo seleccionado Real" << endl;
    }
    else if (solutionType.compare("ArrayReal") == 0)
        solutionType_ = new ArrayRealSolutionType(this) ;
    else {
        cout << "Error: solution type " << solutionType << " invalid" << endl;
        exit(-1) ;
    }

    fx_.resize(numberOfObjectives_);
    x_.resize(numberOfVariables_);
}

ProblemSelf::~ProblemSelf() {
    delete [] lowerLimit_ ;
    delete [] upperLimit_ ;
    delete solutionType_ ;
    std::vector<double>().swap(fx_);
    std::vector<double>().swap(x_);
}

/**
 * Evaluates a solution
 * @param solution The solution to evaluate
 */
void ProblemSelf::evaluate(Solution *solution) {
    std::size_t i;

    XReal * vars = new XReal(solution);

    for (i = 0; i < numberOfVariables_; ++i) {
        x_[i] = vars->getValue(i);
    }

    for (i = 0; i < numberOfObjectives_; ++i){
        fx_[i] = meq_.expb_f(x_)[i];
        solution->setObjective(i, fx_[i]);
    }

    delete vars ;
} // evaluate
