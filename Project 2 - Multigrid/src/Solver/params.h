#pragma once

#include <cmath>

const int VC_V1 = 3;
const int VC_V2 = 3;
const int VC_MAX_ITER = 100;
const double VC_PI = M_PI;
const double VC_EPS = 1e-15;

template <int Dim>
class VC_W{
public:
    static double value;
};

template<>
double VC_W<1>::value = 2./3;

template<>
double VC_W<2>::value = 4./5;
