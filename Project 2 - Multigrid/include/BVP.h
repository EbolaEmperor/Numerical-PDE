#ifndef _BVP_H_
#define _BVP_H_

#include "solver.h"
#include <jsoncpp/json/json.h>
#include <string>

class BVP{
private:
    Json::Value problem;
    rootSolver *solver;

public:
    BVP();
    BVP(const Json::Value & prob);
    BVP(const char * ifile);
    BVP(const std::string & ifile);

    void checkParse();
    void printProblem();
    void read(const char * ifile);
    void read(const std::string & ifile);
    void read(const Json::Value & prob);
    void solve();
    void output(std::ostream & out);
    void output(const std::string & ofile);
    void output();
    std::vector<double> checkError();
};

#endif