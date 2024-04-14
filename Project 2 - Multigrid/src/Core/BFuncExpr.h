#pragma once

#include "Function.h"
#include "mathExpr.h"
#include <string>

template <int Dim>
class BFuncExpr : public Function<Dim>{
private:
    Expression e[2*Dim];

public:
    BFuncExpr();
    BFuncExpr(const BFuncExpr<Dim> &rhs): e(rhs.e) {};
    BFuncExpr(std::shared_ptr<BFuncExpr<Dim>> prhs): BFuncExpr(*prhs) {}

    void setExpr(const std::string &bon, const std::string &expr);

    void addConstant(const std::string &bon, const std::string & name, const double &val);

    double operator () (const double &x) const;

    double operator () (const double &x, const double &y) const;
};