#ifndef _RUNGE_KUTTA_
#define _RUNGE_KUTTA_

#include "IVP.h"

class ClassicalRKSolver : public TimeIntegrator_ConstStep{
protected:
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
};


class ImplicitRKSolver : public TimeIntegrator_ConstStep{
protected:
    virtual ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step) = 0;
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
};


class ESDIRKSolver : public ImplicitRKSolver{
protected:
    ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step);
};


class GaussLegendreRKSolver_OneStepY{
protected:
    std::vector<ColVector> oneStepSolveY(TimeFunction &f, const ColVector &U0, const double &t0, const double &step, const Matrix &A, const ColVector &c);
};


class GaussLegendreRKSolver : public ImplicitRKSolver, public GaussLegendreRKSolver_OneStepY{
protected:
    Matrix A;
    ColVector b, c;
    ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step);
public:
    GaussLegendreRKSolver(const int &stage);
    Matrix getA() const{return A;}
    ColVector getB() const{return b;}
    ColVector getC() const{return c;}
};


class EmbeddedRKSolver : public TimeIntegrator_VariativeStep{
private:
    virtual std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step);
protected:
    int maxorder;
    Matrix A;
    ColVector b1, b2, c;
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps);
};


class FehlbergSolver : public EmbeddedRKSolver{
public:
    FehlbergSolver();
};

class DormandPrinceSolver : public EmbeddedRKSolver{
public:
    DormandPrinceSolver();
};

class AdaptiveGaussLegendreRKSolver : public EmbeddedRKSolver, public GaussLegendreRKSolver_OneStepY{
protected:
    std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step);
public:
    AdaptiveGaussLegendreRKSolver(const int &stage);
};

#endif