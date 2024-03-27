#ifndef _LMM_H_
#define _LMM_H_

#include "IVP.h"

class LMMSolver : public TimeIntegrator_ConstStep{
protected:
    int s, p;
    std::string method;
    virtual ColVector oneStepSolve(TimeFunction &f, const int &i) = 0;
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
};

class AdamsBashforthSolver : public LMMSolver{
private:
    ColVector beta;
protected:
    ColVector oneStepSolve(TimeFunction &f, const int &i);
public:
    AdamsBashforthSolver(const int &order);
};
static void registerAdamsBashforth(void)__attribute__((constructor));


class AdamsMoultonSolver : public LMMSolver{
private:
    ColVector beta;
protected:
    ColVector oneStepSolve(TimeFunction &f, const int &i);
public:
    AdamsMoultonSolver(const int &order);
};
static void registerAdamsMoulton(void)__attribute__((constructor));


class BDFSolver : public LMMSolver{
private:
    ColVector alpha;
    double beta;
protected:
    ColVector oneStepSolve(TimeFunction &f, const int &i);
public:
    BDFSolver(const int &order);
};
static void registerBDF(void)__attribute__((constructor));

#endif