#include "IVP.h"
#include "functionFactory.h"

class HamiltonianFunc : public TimeFunction{
public:
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        const double r = x(0)*x(0) + x(1)*x(1);
        res(0) = x(2);
        res(1) = x(3);
        res(2) = - x(0) / pow(r, 1.5) - 0.015 * x(0) / pow(r, 2.5);
        res(3) = - x(1) / pow(r, 1.5) - 0.015 * x(1) / pow(r, 2.5);
    }
};

static void registerHamiltonian(void)__attribute__((constructor));

void registerHamiltonian(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Hamiltonian", [](double arg){ return (TimeFunction*) new HamiltonianFunc(); });
}

