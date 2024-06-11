#include "IVP.h"
#include "functionFactory.h"

class RestrictedThreeBodyFunc : public TimeFunction{
private:
    double mu;
public:
    RestrictedThreeBodyFunc(const double &_mu): mu(_mu){}
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        res(0) = x(3);
        res(1) = x(4);
        res(2) = x(5);
        res(3) = 2*x(4) + x(0) - mu*(x(0)+mu-1) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*(x(0)+mu) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
        res(4) = -2*x(3) + x(1) - mu*x(1) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*x(1) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
        res(5) = -mu*x(2) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*x(2) / pow(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
    }
    ColVector dp(const ColVector &p, const ColVector &q) const{
        return q;
    }
    ColVector dq(const ColVector &p, const ColVector &q) const{
        ColVector res(3);
        res(0) = 2*q(1) + p(0) - mu*(p(0)+mu-1) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu-1)*(p(0)+mu-1), 1.5) - (1-mu)*(p(0)+mu) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu)*(p(0)+mu), 1.5);
        res(1) = -2*q(0) + p(1) - mu*p(1) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu-1)*(p(0)+mu-1), 1.5) - (1-mu)*p(1) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu)*(p(0)+mu), 1.5);
        res(2) = -mu*p(2) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu-1)*(p(0)+mu-1), 1.5) - (1-mu)*p(2) / pow(p(1)*p(1) + p(2)*p(2) + (p(0)+mu)*(p(0)+mu), 1.5);
        return res;
    }
};

static void registerRestrictedThreeBody(void)__attribute__((constructor));

void registerRestrictedThreeBody(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Restricted 3-Body Problem", [](double arg){ return (TimeFunction*) new RestrictedThreeBodyFunc(arg); });
}
