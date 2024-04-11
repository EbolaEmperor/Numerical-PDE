#pragma once

#include <cmath>

inline double lowBon(const double &x){
    return sin(M_PI*x)/16;
}

struct specialPoint{
    double x, y;
    specialPoint(): x(0.0), y(0.0) {}
    specialPoint(const int &n, const int &i, const int &j){
        double h = 1.0/n;
        double _x = i*h, _y = j*h;
        x = _x;
        y = lowBon(x) + _y * (1.0-lowBon(x));
    }
};
