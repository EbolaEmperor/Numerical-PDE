#pragma once

#include <map>
#include <string>
#include <vector>

class Expression{
private:
    std::string expr;
    std::map<std::string, int> variable;
    std::map<std::string, double> constant;
    int varcnt;

    double evaluate(int &curPos, const std::vector<double> &parm) const;

public:
    Expression():varcnt(0) {}
    Expression(const Expression &rhs):
        expr(rhs.expr), variable(rhs.variable), constant(rhs.constant), varcnt(rhs.varcnt) {}

    void addVariable(const std::string &varName);
    void addConstant(const std::string &conName, const double &val);
    void setExpr(const std::string &s);
    double eval(const std::vector<double> &parm) const;
    double eval() const;
    double operator () (const std::vector<double> &parm) const;
};