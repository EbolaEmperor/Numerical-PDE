#include "IVP.h"
#include "functionFactory.h"

class StabilityFunc : public TimeFunction{
private:
    static double lambda;
public:
    ColVector trueSol(const double &t) const{
        ColVector res(1);
        res(0) = exp(lambda*t);
        return res;
    }
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        res(0) = lambda * x(0);
    }
};

double StabilityFunc::lambda = -1e6;

static void registerStability(void)__attribute__((constructor));

void registerStability(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Stability Problem", [](double arg){ return (TimeFunction*) new StabilityFunc(); });
}
