#pragma once

#include "solver.h"
#include "Boundary.h"

#include "Core/FuncExpr.h"
#include "Core/BFuncExpr.h"

#include <json/json.h>
#include <string>
#include <memory>

template <int Dim>
class BVP{
private:
    Json::Value problem;
    Solver<Dim> *solver;

    std::shared_ptr<Function<Dim>> g;
    FuncExpr<Dim> f;
    BType bonDtil;

public:
    BVP() {};
    BVP(const Json::Value & prob) {read(prob);};
    BVP(const char * ifile) {read(ifile);};
    BVP(const std::string & ifile) {read(ifile);};

    void checkParse();
    void printProblem();
    void read(const char * ifile);
    void read(const std::string & ifile);
    void read(const Json::Value & prob) {problem = prob;}
    void solve();
    void output(std::ostream & out);
    void output(const std::string & ofile);
    void output();
    std::vector<double> checkError();

protected:
    void setFunctions();
};