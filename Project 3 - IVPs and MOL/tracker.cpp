#include <bits/stdc++.h>
#include "CurveTracker.h"
#include "functionFactory.h"
#include "RestrictedThreeBodyFunc.h"
#include "ThreeBodyFunc.h"
#include "StabilityFunc.h"
#include "StiffSampleFunc.h"
#include "VanDerPolFunc.h"
#include "PerturbatedKepler.h"
#include "AliContest.h"
#include "Pendulum.h"
#include <jsoncpp/json/json.h>
using namespace std;

ColVector Circle(double s){
    ColVector res(2);
    res(0) = 0.7 * cos(s * 2 * M_PI);
    res(1) = sin(s * 2 * M_PI) + 2.;
    return res;
}

ColVector Circle2(double s){
    ColVector res(2);
    res(0) = 0.6 * cos(s * 2 * M_PI) + 2 * M_PI;
    res(1) = 0.6 * sin(s * 2 * M_PI) + 1.;
    return res;
}

int main(int argc, char* argv[]){
    Json::Reader reader;
    Json::Value problem;
    ifstream ifs;
    if(argc < 2){
        cerr << "Please provide input file name!" << endl;
        exit(-1);
    }
    ifs.open(argv[1]);
    if(!ifs.is_open()){
        cerr << "cannot read file " << argv[1] << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }

    auto & funcFac = FunctionFactory::Instance();
    TimeFunction *f = funcFac.createFunction(problem["Problem"].asString(), problem["Mass Ratio"].asDouble());

    const double T = problem["T"].asDouble();
    const int order = problem.isMember("Order") ? problem["Order"].asInt() : 0;
    const std::string method = problem["Method"].asString();
    const int nPoints = problem["Points Number"].asInt();
    const int nSteps = problem["Time Steps"].asInt();
    
    CurveTracker tracker;
    tracker.define(method, order, nPoints);
    tracker.solve(*f, Circle, T - 1, nSteps);
    if(problem["Output"].asBool()) tracker.output("result1.txt");

    CurveTracker tracker2;
    tracker2.define(method, order, nPoints);
    tracker2.solve(*f, Circle2, T, nSteps);
    if(problem["Output"].asBool()) tracker2.output("result2.txt");
    return 0;
}
