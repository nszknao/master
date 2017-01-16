#ifndef __EXPFIT_H_INCLUDE_
#define __EXPFIT_H_INCLUDE_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

class MomentEq
{
private:
    std::vector<double> _dG, _lower, _upper;
    std::vector< std::vector<double> > _objWeight;
    std::vector<std::size_t> _objList;
    void _setObjRange();

public:
    std::vector<double> expb_f(const std::vector<double> &);
    static std::vector<double> getMomentFromParameter(const std::vector<double> &);
    void setPrmdG(const std::vector<double> &);
    void setObjList(const std::vector<std::size_t> &);
    void setObjWeight(const std::vector< std::vector<double> > &);
    std::size_t getObjNum();
    std::vector<double> getLowerObj();
    std::vector<double> getUpperObj();
    explicit MomentEq();
    ~MomentEq();
};

#endif // !__EXPFIT_H_INCLUDE_
