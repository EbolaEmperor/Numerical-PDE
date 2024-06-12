#include "CurveTracker.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>

void CurveTracker::define(std::string method, int order, int nPoints){
    tracker.resize(nPoints);
    auto& factory = TimeIntegratorFactory::Instance();
    assert(nPoints > 1);
    for(int i = 0; i < nPoints; ++i)
        tracker[i] = factory.createTimeIntegrator(method, order);
    std::cout << "Curve Tracker with " << nPoints << " points." << std::endl;
}

void CurveTracker::solve(TimeFunction& f, 
                         CurveFunction& curve, 
                         double finalT, int nSteps){
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "Solving... " << std::endl;
    int stTime = clock();
    double h = 1. / (tracker.size() - 1);
    for(int i = 0; i < tracker.size(); ++i){
        ColVector x0 = curve(i * h);
        tracker[i]->solve(f, x0, finalT, nSteps);
    }
    timeSteps = nSteps;
    std::cout << "Solved in " << std::setprecision(3) << (double)(clock()-stTime)/CLOCKS_PER_SEC << "s" << std::endl;
}

void CurveTracker::output(std::string fname){
    std::ofstream fout(fname);
    fout << std::setprecision(12);
    for(int i = 0; i < timeSteps; ++i){
        for(auto& it : tracker){
            auto sol = it->getSolByID(i);
            for(int k = 0; k < sol.size(); ++k)
                fout << sol(k) << " ";
        }
        fout << "\n";
    }
    fout.close();
}