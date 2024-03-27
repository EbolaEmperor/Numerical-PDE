#include "IVP.h"
#include "functionFactory.h"

class HamiltonianFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(4);
        const double r = x(2)*x(2) + x(3)*x(3);
        res(0) = - x(2) / pow(r, 1.5) - 0.015 * x(2) / pow(r, 2.5);
        res(1) = - x(3) / pow(r, 1.5) - 0.015 * x(3) / pow(r, 2.5);
        res(2) = x(0);
        res(3) = x(1);
        return res;
    }
};

static void registerHamiltonian(void)__attribute__((constructor));

void registerHamiltonian(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Hamiltonian", [](double arg){ return (TimeFunction*) new HamiltonianFunc(); });
}