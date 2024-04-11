#pragma once

#include "Core/matrix.h"

template <int Dim>
class IntergridOp{
public:
    virtual ColVector operator () (const ColVector &v) const = 0;
};