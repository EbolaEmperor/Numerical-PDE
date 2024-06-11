#include "IVP.h"
#include "functionFactory.h"

class VanDerPolFunc : public TimeFunction{
public:
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        res(0) = x(1);
        res(1) = 1000*(1-x(0)*x(0))*x(1) - x(0);
    }
};

static void registerVanDerPol(void)__attribute__((constructor));

void registerVanDerPol(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Van der Pol Problem", [](double arg){ return (TimeFunction*) new VanDerPolFunc(); });
}
