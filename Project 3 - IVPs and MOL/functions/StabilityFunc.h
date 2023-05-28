#include "IVP.h"

class StabilityFunc : public TimeFunction{
private:
    static double lambda;
public:
    ColVector trueSol(const double &t) const{
        ColVector res(1);
        res(0) = exp(lambda*t);
        return res;
    }
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(1);
        res(0) = lambda * x(0);
        return res;
    }
};

double StabilityFunc::lambda = -1e6;