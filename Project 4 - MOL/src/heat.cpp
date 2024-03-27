#include <bits/stdc++.h>
#include "IVP.h"
#include <jsoncpp/json/json.h>
using namespace std;

const double PI = acos(-1);

class HeatEquation : public TimeFunction{
private:
    double h;
    double b0(const double &t) const{
        return 0.0;
    }
    double b1(const double &t) const{
        return 0.0;
    }
public:
    HeatEquation(const int &n){
        __isLinear = true;
        h = 1.0 / n;
    }
    ColVector operator () (const ColVector &U, const double &t) const{
        ColVector U1(U.size());
        for(int i = 0; i < U.size(); i++){
            U1(i) = -2*U(i)/(h*h);
            if(i-1 >= 0) U1(i) += U(i-1)/(h*h);
            if(i+1 < U.size()) U1(i) += U(i+1)/(h*h);
            if(i==0) U1(i) += b0(t) / (h*h);
            if(i+1==U.size()) U1(i) += b1(t) / (h*h);
        }
        return U1;
    }
    ColVector solve(const Matrix &a, const ColVector &c, const ColVector &U0, const double &t0, const double &k) const {
        int s = c.size(), m = U0.size();
        ColVector rhs(s*m);
        for(int i = 0; i < s; i++)
            rhs.setSubmatrix(i*m, 0, (*this)(U0, t0+c(i)*k));
        Matrix coef(s*m, s*m);
        for(int i = 0; i < s; i++)
            for(int j = 0; j < s; j++)
                for(int l = 0; l < m; l++){
                    coef(i*m+l, j*m+l) = a(i,j)*k*2.0/(h*h);
                    if(i==j) coef(i*m+l, j*m+l) += 1;
                    if(l) coef(i*m+l, j*m+l-1) = -a(i,j)*k*1.0/(h*h);
                    if(l<m-1) coef(i*m+l, j*m+l+1) = -a(i,j)*k*1.0/(h*h);
                }
        return coef.solve(rhs);
    }
};

double g0(const double &x){
    if(0.45 <= x && x < 0.5) return 20*(x-0.45);
    else if(0.5 <= x && x < 0.55) return -20*(x-0.55);
    else return 0;
}

void outputK(TimeIntegrator *solver, const int &k, const double &step=0){
    char fname[20];
    sprintf(fname, "result%d.txt", k);
    ofstream fout(fname);
    ColVector sol = step ? solver->at(k*step) : solver->getSolByID(k);
    double h = 1.0 / (sol.size()+1);
    fout << "0 0\n";
    for(int i = 0; i < sol.size(); i++){
        fout << (i+1)*h << " " << sol(i) << "\n";
    }
    fout << "1 0\n";
    fout.close();
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

    int n = problem["Space Section"].asInt();
    double h = 1.0 / n;
    double T = problem["End Time"].asDouble();

    ColVector U0(n-1);
    for(int i = 1; i < n; i++)
        U0(i-1) = g0(i*h);
    HeatEquation f(n);

    auto& factory = TimeIntegratorFactory::Instance();
    auto solver = factory.createTimeIntegrator(problem["Method"].asString(), problem["Order"].asInt());
    if(problem.isMember("Time Section")){
        solver->solveWithInfo(f, U0, T, problem["Time Section"].asInt());
    } else {
        if(problem.isMember("Max Step"))
            solver->setMaxStep(problem["Max Step"].asDouble());
        solver->solveWithInfo(f, U0, T, problem["Tolerance"].asDouble());
    }
    
    if(problem["Output"].asBool()){
        // Example中的要求
        double step = problem.isMember("Time Section") ? 0 : T/10;
        outputK(solver, 1, step);
        outputK(solver, 2, step);
        outputK(solver, 10, step);
    }
    if(problem["Dense-Discrete Output"].asBool())
        solver->denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble());
    return 0;
}