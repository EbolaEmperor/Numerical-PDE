#pragma once

#include "Function.h"
#include "mathExpr.h"

template <int Dim>
class FuncExpr : public Function<Dim>{
private:
    Expression e;
    Expression Delta_e;

public:
    FuncExpr();
    FuncExpr(const FuncExpr<Dim>& rhs): e(rhs.e), Delta_e(rhs.e) {}
    FuncExpr(std::shared_ptr<FuncExpr<Dim>> prhs): FuncExpr(*prhs) {}

    void setExpr(const std::string &expr){
        e.setExpr(expr);
    }

    double operator () (const double &x) const;

    void setDeltaExpr(const std::string &expr);

    double operator () (const double &x, const double &y) const;
    
    double delta(const double &x, const double &y) const;
};