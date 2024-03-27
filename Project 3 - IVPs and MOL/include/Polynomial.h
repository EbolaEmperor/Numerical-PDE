#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "matrix.h"

// 一个不完整的多项式运算库

class Polynomial{
private:
    int n;
    double *a;
    void removeLeadingZero();
public:
    Polynomial(): n(-1), a(nullptr){}
    Polynomial(const int &_n);
    Polynomial(const int &_n, const double *p);
    Polynomial(const Polynomial &rhs);
    Polynomial& operator = (const Polynomial &rhs);
    Polynomial operator + (const Polynomial &rhs) const;
    Polynomial operator * (const Polynomial &rhs) const;
    friend Polynomial pow(Polynomial lhs, int b);
    Polynomial operator / (const double &c) const;
    Polynomial operator / (const Polynomial &rhs) const;
    Polynomial derivative() const;
    Polynomial integral() const;
    std::vector<Complex> roots();
    double operator () (const double &x) const;
    double& coef(const int &i);
    double coef(const int &i) const;
    void print() const;
    friend std::ostream& operator << (std::ostream &out, const Polynomial &p);
};

Polynomial constPolynomial(const double &a);
Polynomial linearPolynomial(const double &a, const double &b);
Polynomial HermiteInterpolation32(const double &x0, const double &x1, const double &y0, const double &dy0, const double &y1, const double &dy1);

#endif