#include "IVP.h"

class VanDerPolFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(2);
        res(0) = x(1);
        res(1) = 1000*(1-x(0)*x(0))*x(1) - x(0);
        return res;
    }
};