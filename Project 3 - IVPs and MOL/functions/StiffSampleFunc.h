#include "IVP.h"
#include "functionFactory.h"

class StiffSampleFunc : public TimeFunction{
private:
    static double lambda, eta;
public:
    ColVector trueSol (const ColVector &x, const double &t) const{
        ColVector res(1);
        res(0) = exp(lambda * t) * (eta - 1) + cos(t);
        return res;
    }
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        res(0) = lambda * ( x(0) - cos(t) ) - sin(t);
    }
};

double StiffSampleFunc::lambda = -1e6;
double StiffSampleFunc::eta = 1.5;

static void registerStiffSample(void)__attribute__((constructor));

void registerStiffSample(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Stiff Sample Problem", [](double arg){ return (TimeFunction*) new StiffSampleFunc(); });
}
