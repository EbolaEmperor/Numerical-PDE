#include <bits/stdc++.h>
#include "matrix.h"
#include <jsoncpp/json/json.h>
using namespace std;

double mu;

ColVector FTCS(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    for(int i = 0; i < m; i++){
        U1(i) = U0(i);
        U1(i) += mu/2 * U0(i ? i-1 : m-1);
        U1(i) -= mu/2 * U0(i<m-1 ? i+1 : 0);
    }
    return U1;
}

ColVector leapfrog(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    static ColVector Upre;
    if(Upre.size()==0) Upre = U0;
    for(int i = 0; i < m; i++){
        U1(i) = Upre(i);
        U1(i) += mu * U0(i ? i-1 : m-1);
        U1(i) -= mu * U0(i<m-1 ? i+1 : 0);
    }
    Upre = U0;
    return U1;
}

ColVector LaxFriedrichs(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    for(int i = 0; i < m; i++){
        U1(i) += (1+mu)/2 * U0(i ? i-1 : m-1);
        U1(i) += (1-mu)/2 * U0(i<m-1 ? i+1 : 0);
    }
    return U1;
}

ColVector LaxWendroff(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    for(int i = 0; i < m; i++){
        U1(i) = (1-mu*mu) * U0(i);
        U1(i) += mu*(1+mu)/2 * U0(i ? i-1 : m-1);
        U1(i) += mu*(-1+mu)/2 * U0(i<m-1 ? i+1 : 0);
    }
    return U1;
}

ColVector upwind(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    for(int i = 0; i < m; i++){
        if(mu>=0){
            U1(i) = (1-mu) * U0(i);
            U1(i) += mu * U0(i ? i-1 : m-1);
        } else {
            U1(i) = (1+mu) * U0(i);
            U1(i) -= mu * U0(i<m-1 ? i+1 : 0);
        }
    }
    return U1;
}

ColVector BeamWarming(const ColVector &U0){
    int m = U0.size();
    ColVector U1(m);
    for(int i = 0; i < m; i++){
        if(mu>=0){
            U1(i) = (1 - 1.5*mu + 0.5*mu*mu) * U0(i);
            U1(i) += (2*mu - mu*mu) * U0(i ? i-1 : m-1);
            U1(i) += (-0.5*mu + 0.5*mu*mu) * U0(i-1>0 ? i-2 : m-2);
        } else {
            U1(i) = (1 + 1.5*mu + 0.5*mu*mu) * U0(i);
            U1(i) += (-2*mu - mu*mu) * U0(i<m-1 ? i+1 : 0);
            U1(i) += (0.5*mu + 0.5*mu*mu) * U0(i<m-2 ? i+2 : 1);
        }
    }
    return U1;
}

map<string, ColVector (*) (const ColVector &)> oneStepMethods;

void registerMethods(){
    oneStepMethods["FTCS"] = FTCS;
    oneStepMethods["leapfrog"] = leapfrog;
    oneStepMethods["Lax-Friedrichs"] = LaxFriedrichs;
    oneStepMethods["Lax-Wendroff"] = LaxWendroff;
    oneStepMethods["upwind"] = upwind;
    oneStepMethods["Beam-Warming"] = BeamWarming;
}

// Initial condition
double f0(const double &x){
    return exp(-20*sqr(x-2)) + exp(-sqr(x-5));
}

vector<ColVector> solveAdvection(const string &method, const double &l, const double &r, const double &h, const double &T, const double &k){
    int m = (int)round((r-l)/h);
    ColVector U(m+1);
    vector<ColVector> sol;
    sol.push_back(U);
    for(int i = 0; i <= m; i++){
        U(i) = f0(l + i*h);
    }
    if(!oneStepMethods.count(method)){
        cerr << "No such method called '" << method << "'" << endl;
        exit(-1);
    }
    mu = k / h; // We fix a=1
    auto fun = oneStepMethods[method];
    for(double t = 0; t < T - 1e-6; t += k){
        U = fun(U);
        sol.push_back(U);
    }
    return sol;
}

int main(int argc, char* argv[]){
    registerMethods();
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
    
    cout << "Solving with " << problem["Method"].asString() << " method..." << endl;
    auto res = solveAdvection(problem["Method"].asString(),
                              problem["Left"].asDouble(),
                              problem["Right"].asDouble(),
                              problem["Space Step"].asDouble(),
                              problem["End Time"].asDouble(),
                              problem["Time Step"].asDouble());
    cout << "Solved." << endl;

    ofstream fout("result.txt");
    double k = problem["Time Step"].asDouble();
    int skip = problem.isMember("Output Skip") ? problem["Output Skip"].asInt() : 1;
    for(int i = 0; i < res.size(); i += skip){
        fout << i*k << " ";
        for(int j = 0; j < res[i].size(); j++)
            fout << res[i](j) << " ";
        fout << "\n";
    }
    fout.close();
    cout << "Result has been saved to result.txt." << endl;

    return 0;
}