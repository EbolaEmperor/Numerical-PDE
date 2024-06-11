#include "IVP.h"
#include "functionFactory.h"

class AliContestFunc : public TimeFunction{
private:
    double omega;
public:
    AliContestFunc(const double &arg){
        omega = arg;
    }
    void compute (const ColVector &x, const double &t, ColVector &res) const{
        res(0) = x(1);
        res(1) = -(1 + 0.01 * cos(omega*t)) * x(0);
    }
};

static void registerAliContest(void)__attribute__((constructor));

void registerAliContest(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Alibaba Contest", [](double arg){ return (TimeFunction*) new AliContestFunc(arg); });
}
