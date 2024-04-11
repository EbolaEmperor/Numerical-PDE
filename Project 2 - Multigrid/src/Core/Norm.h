#pragma once

#include <iostream>
#include <cmath>
#include <vector>

class Norm{
public:
    virtual double operator () (const std::vector<double> &vec) const = 0;
};

class Norm_p : public Norm{
private:
    double p;
public:
    Norm_p(const int &p): p(p) {}
    double operator () (const std::vector<double> &vec) const{
        double sum = 0;
        for(const double &x : vec)
            sum += std::pow(fabs(x), p);
        sum /= vec.size();
        return pow(sum, 1.0/p);
    }
};

class Norm_inf : public Norm{
public:
    double operator () (const std::vector<double> &vec) const{
        double mx = 0;
        for(const double &x : vec)
            mx = std::max(mx, fabs(x));
        return mx;
    }
};