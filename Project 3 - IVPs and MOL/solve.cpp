#include <bits/stdc++.h>
#include "IVP.h"
#include "RestrictedThreeBodyFunc.h"
#include "ThreeBodyFunc.h"
#include "StabilityFunc.h"
#include "StiffSampleFunc.h"
#include "VanDerPolFunc.h"
#include "Custom.h"
#include "json.h"
using namespace std;

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

    TimeFunction *f;
    if(problem["Problem"].asString() == "Restricted 3-Body Problem"){
        f = new RestrictedThreeBodyFunc(problem["Mass Ratio"].asDouble());
    } else if(problem["Problem"].asString() == "3-Body Problem"){
        f = new ThreeBodyFunc();
    } else if(problem["Problem"].asString() == "Van der Pol Problem"){
        f = new VanDerPolFunc();
    } else if(problem["Problem"].asString() == "Stability Problem"){
        f = new StabilityFunc();
    } else if(problem["Problem"].asString() == "Stiff Sample Problem"){
        f = new StiffSampleFunc();
    } else if(problem["Problem"].asString() == "Custom"){
        f = new CustomFunc();
    } else {
        cerr << "[Error] No such problem called " << problem["Problem"].asString() << endl;
        exit(-1);
    }

    const int m = problem["Init"].size();
    ColVector x0(m);
    for(int i = 0; i < m; i++)
        x0(i) = problem["Init"][i].asDouble();
    const double T = problem["T"].asDouble();
    const int periodic = problem.isMember("Periodic") ? problem["Periodic"].asInt() : 1;
    const int order = problem.isMember("Order") ? problem["Order"].asInt() : 0;

    auto& factory = TimeIntegratorFactory::Instance();
    auto solver = factory.createTimeIntegrator(problem["Method"].asString(), order);
    TimeIntegrator* solver2 = nullptr;

    if(problem.isMember("Max Step")){
        solver->setMaxStep(problem["Max Step"].asDouble());
    }
    if(problem.isMember("eps"))
        solver->solveWithInfo(*f, x0, T*periodic, problem["eps"].asDouble());
    else
        solver->solveWithInfo(*f, x0, T*periodic, problem["Grid Size"].asInt());
    
    if(problem["Output"].asBool())
        solver->output("result.txt");
    if(problem["Dense Output"].asBool())
        solver->denseOutput("result-dense.txt");
    if(problem["Dense-Discrete Output"].asBool()){
        const double step = problem.isMember("Dense-Discrete Output Step") ? problem["Dense-Discrete Output Step"].asDouble() : 0.01;
        solver->denseDiscreteOutput("result-dense-discrete.txt", T/ceil(T/step));
    }

    if( (problem["Convergence Analysis"].asBool() || problem["Richardson Error Estimate"].asBool()) && problem.isMember("Grid Size") ){
        cout << "Resolving with " << (problem["Grid Size"].asInt()<<1) << " steps for error-estimation or convergence-analysis..." << endl;
        solver2 = factory.createTimeIntegrator(problem["Method"].asString(), order);
        solver2->solve(*f, x0, T*periodic, problem["Grid Size"].asInt()<<1);
    }

    if(problem["Richardson Error Estimate"].asBool() && problem.isMember("Grid Size")){
        cout << "--------------------------------------------------------------------------------" << endl;
        cout << "Error estimate with Richardson extrapolation." << endl;
        auto err1 = abs(solver->at(periodic*T)-solver2->at(periodic*T));
        cout << "Error Estimate: " << err1.T() << endl;
        if(problem["Convergence Analysis"].asBool()){
            cout << "Resolving with " << (problem["Grid Size"].asInt()<<2) << " steps for convergence-analysis..." << endl;
            auto solver4 = factory.createTimeIntegrator(problem["Method"].asString(), order);
            solver4->solve(*f, x0, T*periodic, problem["Grid Size"].asInt()<<2);
            err1 = abs(solver->at(periodic*T)-solver4->at(periodic*T));
            auto err2 = abs(solver2->at(periodic*T)-solver4->at(periodic*T));
            cout << "Convergence rate: " << log2(dotdiv(err1,err2)-1.0).T() << endl;
        }
    }
    
    if(problem.isMember("Periodic")){
        cout << "--------------------------------------------------------------------------------" << endl;
        auto err1 = abs(solver->at(periodic*T)-x0);
        cout << "Error between 0 and t_end: " << err1.T() << endl;
        if(problem["Convergence Analysis"].asBool() && problem.isMember("Grid Size")){
            auto err2 = abs(solver2->at(periodic*T)-x0);
            cout << "Convergence rate: " << log2(dotdiv(err1,err2)).T() << endl;
        }
    }
    return 0;
}
