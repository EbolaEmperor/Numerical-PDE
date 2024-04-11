#include "BFuncExpr.h"

#include "Irregular.h"

#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
using std::cerr;
using std::endl;
using std::vector;

const double BFE_eps = 1e-15;

template<>
BFuncExpr<1>::BFuncExpr(){}

template<>
void BFuncExpr<1>::setExpr(const std::string &bon, const std::string &expr){
    if(bon == "x=0") e[0].setExpr(expr);
    else if(bon == "x=1") e[1].setExpr(expr);
    else {
        cerr << "[Error] Unrecognized bondary '" << bon << "'" << endl;
        exit(-1);
    }
}

template<>
double BFuncExpr<1>::operator () (const double &x) const{
    vector<double> parm;
    if(fabs(x) < BFE_eps) return e[0](parm);
    else if(fabs(x-1) < BFE_eps) return e[1](parm);
    else{
        cerr << "[Error] Point (" << x << ") is not on the bondary" << endl;
        exit(-1);
    } 
}

template<>
double BFuncExpr<2>::operator () (const double &x) const{
    assert(false); return 0.;
}

template<>
BFuncExpr<2>::BFuncExpr(){
    for(int i = 0; i < 4; i++){
        e[i].addVariable("x");
        e[i].addVariable("y");
    }
}

template<>
void BFuncExpr<2>::addConstant(const std::string &bon, const std::string & name, const double &val){
    if(bon == "x=0") e[0].addConstant(name, val);
    else if(bon == "x=1") e[1].addConstant(name, val);
    else if(bon == "Down Bondary") e[2].addConstant(name, val);
    else if(bon == "y=1") e[3].addConstant(name, val);
    else {
        cerr << "[Error] Unrecognized bondary '" << bon << "'" << endl;
        exit(-1);
    }
}

template<>
void BFuncExpr<2>::setExpr(const std::string &bon, const std::string &expr){
    if(bon == "x=0") e[0].setExpr(expr);
    else if(bon == "x=1") e[1].setExpr(expr);
    else if(bon == "Down Bondary") e[2].setExpr(expr);
    else if(bon == "y=1") e[3].setExpr(expr);
    else {
        cerr << "[Error] Unrecognized bondary '" << bon << "'" << endl;
        exit(-1);
    }
}

template<>
double BFuncExpr<1>::operator () (const double &x, const double &y) const{
    assert(false); return 0.;
}

template<>
double BFuncExpr<2>::operator () (const double &x, const double &y) const{
    vector<double> parm;
    parm.push_back(x);
    parm.push_back(y);
    if(fabs(y) < BFE_eps) return e[2](parm);
    else if(fabs(y-1) < BFE_eps) return e[3](parm);
    else if(fabs(x) < BFE_eps) return e[0](parm);
    else if(fabs(x-1) < BFE_eps) return e[1](parm);
    else if(fabs(y-lowBon(x)) < BFE_eps) return e[2](parm);
    else{
        cerr << "[Error] Point (" << x << "," << y << ") is not on the bondary" << endl;
        exit(-1);
    } 
}

template class BFuncExpr<1>;
template class BFuncExpr<2>;