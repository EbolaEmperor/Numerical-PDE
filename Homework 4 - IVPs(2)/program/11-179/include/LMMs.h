#ifndef _LMM_H_
#define _LMM_H_

#include "IVP.h"

class AdamsBashforthSolver : public TimeIntegrator_ConstStep{
private:
    ColVector beta;
protected:
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
public:
    AdamsBashforthSolver(const int &order);
};

class AdamsMoultonSolver : public TimeIntegrator_ConstStep{
private:
    ColVector beta;
protected:
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
public:
    AdamsMoultonSolver(const int &order);
};

class BDFSolver : public TimeIntegrator_ConstStep{
private:
    ColVector alpha;
    double beta;
protected:
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
public:
    BDFSolver(const int &order);
};

#endif