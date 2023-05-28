#include "IVP.h"

class StiffSampleFunc : public TimeFunction{
private:
    static double lambda, eta;
public:
    ColVector trueSol (const ColVector &x, const double &t) const{
        ColVector res(1);
        res(0) = exp(lambda * t) * (eta - 1) + cos(t);
        return res;
    }
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(1);
        res(0) = lambda * ( x(0) - cos(t) ) - sin(t);
        return res;
    }
};

double StiffSampleFunc::lambda = -1e6;
double StiffSampleFunc::eta = 1.5;