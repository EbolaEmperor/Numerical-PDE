#pragma once

#include "Core/matrix.h"
#include "IntergridOp.h"

template <int Dim>
class FullWeighting : public IntergridOp<Dim>{
public:
    ColVector operator () (const ColVector &v) const;
};