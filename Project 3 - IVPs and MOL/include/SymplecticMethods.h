#ifndef _SYMPLECTIC_METHODS_H_
#define _SYMPLECTIC_METHODS_H_

#include "IVP.h"

class SymplecticSolver : public TimeIntegrator_ConstStep{
protected:
    virtual ColVector oneStepSolve(TimeFunction &f, 
                                   const ColVector &x0, 
                                   const ColVector &f0, 
                                   const double &t0, 
                                   const double &step) = 0;
public:
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
};


class SymplecticEulerSolver : public SymplecticSolver{
protected:
    virtual ColVector oneStepSolve(TimeFunction &f, 
                                   const ColVector &x0, 
                                   const ColVector &f0, 
                                   const double &t0, 
                                   const double &step);
public:
    void oneStepP(TimeFunction &f, 
                  const ColVector &p0,
                  const ColVector &q0,
                  ColVector &p1,
                  ColVector &q1,
                  const double &step);
    void oneStepQ(TimeFunction &f, 
                  const ColVector &p0,
                  const ColVector &q0,
                  ColVector &p1,
                  ColVector &q1,
                  const double &step);
};
static void registerSymplecticEuler(void)__attribute__((constructor));


class StormerVerletSolver : public SymplecticSolver{
protected:
    SymplecticEulerSolver symEuler;
    virtual ColVector oneStepSolve(TimeFunction &f, 
                                   const ColVector &x0, 
                                   const ColVector &f0, 
                                   const double &t0, 
                                   const double &step);
};
static void registerStormerVerlet(void)__attribute__((constructor));

#endif