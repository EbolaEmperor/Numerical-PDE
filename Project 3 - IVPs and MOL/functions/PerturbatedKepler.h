#include "IVP.h"
#include "functionFactory.h"

class PerturbatedKeplerFunc : public TimeFunction{
public:
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        const double r = x(0)*x(0) + x(1)*x(1);
        res(0) = x(2);
        res(1) = x(3);
        res(2) = - x(0) / pow(r, 1.5) - 0.015 * x(0) / pow(r, 2.5);
        res(3) = - x(1) / pow(r, 1.5) - 0.015 * x(1) / pow(r, 2.5);
    }
    ColVector dp(const ColVector &p, const ColVector &q) const{
        return q;
    }
    ColVector dq(const ColVector &p, const ColVector &q) const{
        const double r = p(0)*p(0) + p(1)*p(1);
        ColVector res(2);
        res(0) = - p(0) / pow(r, 1.5) - 0.015 * p(0) / pow(r, 2.5);
        res(1) = - p(1) / pow(r, 1.5) - 0.015 * p(1) / pow(r, 2.5);
        return res;
    }
};

static void registerPerturbatedKepler(void)__attribute__((constructor));

void registerPerturbatedKepler(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Perturbated Kepler", [](double arg){ return (TimeFunction*) new PerturbatedKeplerFunc(); });
}

