#include "IVP.h"
#include "functionFactory.h"

class PendulumFunc : public TimeFunction{
public:
    void compute(const ColVector &x, const double &t, ColVector &res) const{
        res(0) = x(1);
        res(1) = -sin(x(0));
    }
};

static void registerPendulum(void)__attribute__((constructor));

void registerPendulum(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Pendulum", [](double arg){ return (TimeFunction*) new PendulumFunc(); });
}
