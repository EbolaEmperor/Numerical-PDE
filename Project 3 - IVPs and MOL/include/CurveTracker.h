#ifndef _CURVE_TRACKER_H_
#define _CURVE_TRACKER_H_

#include "IVP.h"
#include <vector>
#include <string>
#include <functional>

typedef ColVector CurveFunction(double s);

class CurveTracker{
protected:
    int timeSteps;
    std::vector<TimeIntegrator*> tracker;
public:
    void define(std::string method, int order, int nPoints);
    void solve(TimeFunction& f, CurveFunction& curve, double finalT, int nSteps);
    void output(std::string fname);
};

#endif