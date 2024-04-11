#pragma once

template <int Dim>
class Function{
public:
    virtual double operator () (const double &x) const {return 0.;};
    virtual double operator () (const double &x, const double &y) const {return 0.;};
    virtual double delta(const double &x, const double &y) const {return 0.;};
};