#include "IVP.h"
#include "functionFactory.h"

class ThreeBodyFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(18);
        res(0) = x(9);
        res(1) = x(10);
        res(2) = x(11);
        res(3) = x(12);
        res(4) = x(13);
        res(5) = x(14);
        res(6) = x(15);
        res(7) = x(16);
        res(8) = x(17);
        res(9) = (x(3)-x(0))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(6)-x(0))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5);
        res(10) = (x(4)-x(1))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(7)-x(1))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5);
        res(11) = (x(5)-x(2))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(8)-x(2))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5);
        res(12) = (x(0)-x(3))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(6)-x(3))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        res(13) = (x(1)-x(4))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(7)-x(4))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        res(14) = (x(2)-x(5))/pow( sqr(x(3)-x(0)) + sqr(x(4)-x(1)) + sqr(x(5)-x(2)), 1.5) + (x(8)-x(5))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        res(15) = (x(0)-x(6))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5) + (x(3)-x(6))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        res(16) = (x(1)-x(7))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5) + (x(4)-x(7))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        res(17) = (x(2)-x(8))/pow( sqr(x(6)-x(0)) + sqr(x(7)-x(1)) + sqr(x(8)-x(2)), 1.5) + (x(5)-x(8))/pow( sqr(x(6)-x(3)) + sqr(x(7)-x(4)) + sqr(x(8)-x(5)), 1.5);
        return res;
    }
};

static void registerThreeBody(void)__attribute__((constructor));

void registerThreeBody(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("3-Body Problem", [](double arg){ return (TimeFunction*) new ThreeBodyFunc(); });
}