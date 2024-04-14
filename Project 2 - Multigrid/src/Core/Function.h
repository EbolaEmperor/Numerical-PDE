#pragma once

#include <memory>

template <int Dim>
class Function{
public:
    Function() {}
    Function(const Function<Dim>&) {}
    Function(std::shared_ptr<Function<Dim>> prhs) {}

    virtual double operator () (const double &x) const {return 0.;};
    virtual double operator () (const double &x, const double &y) const {return 0.;};
    virtual double delta(const double &x, const double &y) const {return 0.;};
};