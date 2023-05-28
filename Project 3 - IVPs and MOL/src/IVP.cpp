#include "IVP.h"
#include "Polynomial.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <iomanip>
using namespace std;

//--------------------TimeIntegrator Factory----------------------------

bool TimeIntegratorFactory::registerTimeIntegrator(const string &ID, CreateTimeIntegratorCallBack createFn){
    if(callbacks_.count(ID)) return false;
    callbacks_[ID] = createFn;
    return true;
}

bool TimeIntegratorFactory::unregisterTimeIntegrator(const string &ID){
    if(!callbacks_.count(ID)) return false;
    callbacks_.erase(ID);
    return true;
}

TimeIntegrator* TimeIntegratorFactory::createTimeIntegrator(const string &ID, const int &order){
    if(!callbacks_.count(ID)){
        cerr << "TimeIntegratorFactory:: No such TimeIntegrator called '" << ID << "'." << endl;
        return nullptr;
    }
    return callbacks_[ID](order);
}

TimeIntegratorFactory& TimeIntegratorFactory::Instance(){
    static TimeIntegratorFactory factory;
    return factory;
}

//-------------------TimeIntegrator_ConstStep----------------------------

ColVector TimeIntegrator_ConstStep::at(const double &t) const{
    if(t < 0 || t > maxTime + 1){
        cerr << "[Error] TimeIntegrator at: out of range." << endl;
        exit(-1);
    }
    int i = floor(t/timeStep);
    if(i==sol.size()-1) return sol[i] + (t-i*timeStep) * dsol[i];
    int m = sol[0].size();
    ColVector res(m);
    for(int j = 0; j < m; j++){
        Polynomial p = HermiteInterpolation32(i*timeStep, (i+1)*timeStep, sol[i](j), dsol[i](j), sol[i+1](j), dsol[i+1](j));
        res(j) = p(t);
    }
    return res;
}

void TimeIntegrator_ConstStep::output(const string &fname) const{
    cout << "--------------------------------------------------------------------------------" << endl;
    ofstream fout(fname);
    fout << setprecision(12);
    for(int i = 0; i < sol.size(); i++){
        fout << i*timeStep;
        for(int j = 0; j < sol[i].size(); j++)
            fout << ' ' << sol[i](j);
        fout << '\n';
    }
    fout.close();
    cout << "Point Output: Results has been saved to " << fname << endl;
}

void TimeIntegrator_ConstStep::denseOutput(const string &fname) const{
    cout << "--------------------------------------------------------------------------------" << endl;
    ofstream fout(fname);
    fout << setprecision(12);
    const int m = sol[0].size();
    for(int i = 0; i < sol.size()-1; i++){
        fout << i*timeStep;
        for(int j = 0; j < m; j++){
            Polynomial p = HermiteInterpolation32(i*timeStep, (i+1)*timeStep, sol[i](j), dsol[i](j), sol[i+1](j), dsol[i+1](j));
            fout << " " << p;
        }
        fout << '\n';
    }
    fout.close();
    cout << "Dense Output: Results has been saved to " << fname << endl;
}

void TimeIntegrator_ConstStep::solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps){
    cerr << "[Error] You should provide GridSize at a const-step solver." << endl;
    exit(-1);
}

//--------------------------TimeIntegrator_VariativeStep---------------------------------

ColVector TimeIntegrator_VariativeStep::at(const double &t) const{
    if(t < 0 || t > maxTime + 1){
        cerr << "[Error] TimeIntegrator at: out of range." << endl;
        exit(-1);
    }
    int i = upper_bound(timePoint.begin(), timePoint.end(), t) - timePoint.begin() - 1;
    if(i==sol.size()-1) return sol[i] + (t-timePoint[i]) * dsol[i];
    int m = sol[0].size();
    ColVector res(m);
    for(int j = 0; j < m; j++){
        Polynomial p = HermiteInterpolation32(timePoint[i], timePoint[i+1], sol[i](j), dsol[i](j), sol[i+1](j), dsol[i+1](j));
        res(j) = p(t);
    }
    return res;
}

void TimeIntegrator_VariativeStep::output(const string &fname) const{
    cout << "--------------------------------------------------------------------------------" << endl;
    ofstream fout(fname);
    fout << setprecision(12);
    for(int i = 0; i < sol.size(); i++){
        fout << timePoint[i];
        for(int j = 0; j < sol[i].size(); j++)
            fout << ' ' << sol[i](j);
        fout << '\n';
    }
    fout.close();
    cout << "Point Output: Results has been saved to " << fname << endl;
}

void TimeIntegrator_VariativeStep::denseOutput(const string &fname) const{
    cout << "--------------------------------------------------------------------------------" << endl;
    ofstream fout(fname);
    fout << setprecision(12);
    const int m = sol[0].size();
    for(int i = 0; i < sol.size()-1; i++){
        fout << timePoint[i];
        for(int j = 0; j < m; j++){
            Polynomial p = HermiteInterpolation32(timePoint[i], timePoint[i+1], sol[i](j), dsol[i](j), sol[i+1](j), dsol[i+1](j));
            fout << ' ' << p;
        }
        fout << '\n';
    }
    fout.close();
    cout << "Dense Output: Results has been saved to " << fname << endl;
}

void TimeIntegrator_VariativeStep::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize){
    cerr << "[Error] You cannot provide GridSize at a variative-step solver." << endl;
    exit(-1);
}

//---------------------------------Dense Discrete Output-----------------------------------

void TimeIntegrator::denseDiscreteOutput(const std::string &fname, const double &step) const{
    cout << "--------------------------------------------------------------------------------" << endl;
    ofstream fout(fname);
    fout << setprecision(12);
    for(double cur = 0; cur <= maxTime + step/2; cur += step){
        fout << cur;
        ColVector y = at(cur);
        for(int j = 0; j < y.size(); j++)
            fout << ' ' << y(j);
        fout << '\n';
    }
    fout.close();
    cout << "Dense Discrete Output: Results has been saved to " << fname << endl;
}

void TimeIntegrator::solveWithInfo(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize){
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Solving... " << endl;
    int stTime = clock();
    solve(f, x0, T, gridSize);
    cout << "Solved in " << setprecision(3) << (double)(clock()-stTime)/CLOCKS_PER_SEC << "s" << endl;
    cout << setprecision(6) << "Solution at t_end: " << sol.back().T() << endl;
}

void TimeIntegrator::solveWithInfo(TimeFunction &f, const ColVector &x0, const double &T, const double &eps){
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Solving... (with error tolerance " << eps << ")" << endl;
    int stTime = clock();
    solve(f, x0, T, eps);
    cout << "Solved in " << setprecision(3) << (double)(clock()-stTime)/CLOCKS_PER_SEC << "s" << endl;
    cout << "Grid Size: " << sol.size()-1 << endl;
    cout << setprecision(6) << "Solution at t_end: " << sol.back().T() << endl;
}

void TimeIntegrator::setMaxStep(const double &t){
    maxStep = t;
}

ColVector TimeIntegrator::getSolByID(const int &i) const{
    if(i < 0 || i >= sol.size()){
        std::cerr << "[Error] getSolByID:: Out of range!" << std::endl;
        exit(-1);
    }
    return sol[i];
}