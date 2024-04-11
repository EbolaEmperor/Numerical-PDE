#include "FuncExpr.h"

#include <cassert>
#include <vector>
using std::vector;

template<>
FuncExpr<1>::FuncExpr(){
    e.addVariable("x");
}

template<>
FuncExpr<2>::FuncExpr(){
    e.addVariable("x");
    e.addVariable("y");
}

template<>
double FuncExpr<1>::operator () (const double &x) const{
    vector<double> parm;
    parm.push_back(x);
    return e(parm);
}

template<>
void FuncExpr<1>::setDeltaExpr(const std::string &expr){
    assert(false);
}

template<>
double FuncExpr<1>::operator () (const double &x, const double &y) const{
    assert(false); return 0.;
}

template<>
double FuncExpr<1>::delta(const double &x, const double &y) const{
    assert(false); return 0.;
}

template<>
double FuncExpr<2>::operator () (const double &x) const{
    assert(false); return 0.;
}


template<>
void FuncExpr<2>::setDeltaExpr(const std::string &expr){
    Delta_e.addVariable("x");
    Delta_e.addVariable("y");
    Delta_e.setExpr(expr);
}

template<>
double FuncExpr<2>::operator () (const double &x, const double &y) const{
    vector<double> parm;
    parm.push_back(x);
    parm.push_back(y);
    return e(parm);
}

template<>
double FuncExpr<2>::delta(const double &x, const double &y) const{
    vector<double> parm;
    parm.push_back(x);
    parm.push_back(y);
    return Delta_e(parm);
}

template class FuncExpr<1>;
template class FuncExpr<2>;