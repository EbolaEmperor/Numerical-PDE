#include "SymplecticMethods.h"

//-----------------------------------Base Tools-----------------------------------------

ColVector hamiltonianGetP(const ColVector &x){
    ColVector res(x.size() / 2);
    for(int i = 0; i < x.size() / 2; ++i)
        res(i) = x(i);
    return res;
}

ColVector hamiltonianGetQ(const ColVector &x){
    ColVector res(x.size() / 2);
    for(int i = 0; i < x.size() / 2; ++i)
        res(i) = x(i + x.size() / 2);
    return res;
}

ColVector hamiltonianMerge(const ColVector &p, const ColVector &q){
    ColVector res(p.size() + q.size());
    for(int i = 0; i < p.size(); ++i) res(i) = p(i);
    for(int i = 0; i < q.size(); ++i) res(i + p.size()) = q(i);
    return res;
}

//-----------------------------------Symplectic Solver-----------------------------------------

void SymplecticSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], dsol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i], i * timeStep));
    }
}

//-----------------------------------Symplectic Euler Method-----------------------------------------

void registerSymplecticEuler(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Symplectic Euler", [](int p){ return (TimeIntegrator*) new SymplecticEulerSolver(); });
}

ColVector SymplecticEulerSolver::oneStepSolve(TimeFunction &f, 
                                   const ColVector &x0, 
                                   const ColVector &f0, 
                                   const double &t0, 
                                   const double &step){
    auto p0 = hamiltonianGetP(x0);
    auto q0 = hamiltonianGetQ(x0);
    auto p1 = p0; auto q1 = q0;
    oneStepP(f, p0, q0, p1, q1, step);
    return hamiltonianMerge(p1, q1);
}

void SymplecticEulerSolver::oneStepP(TimeFunction &f, 
                  const ColVector &p0,
                  const ColVector &q0,
                  ColVector &p1,
                  ColVector &q1,
                  const double &step){
    auto ptmp = p0; p1 = p0;
    int cnt = 0;
    do{
        ptmp = p1;
        p1 = p0 + step * f.dp(p1, q0);
    } while(++cnt < 100 && relativeError(ptmp, p1) > 1e-14);
    q1 = q0 + step * f.dq(p1, q0);
}

void SymplecticEulerSolver::oneStepQ(TimeFunction &f, 
                  const ColVector &p0,
                  const ColVector &q0,
                  ColVector &p1,
                  ColVector &q1,
                  const double &step){
    auto qtmp = q0; q1 = q0;
    int cnt = 0;
    do{
        qtmp = q1;
        q1 = q0 + step * f.dq(p0, q1);
    } while(++cnt < 100 && relativeError(qtmp, q1) > 1e-14);
    p1 = p0 + step * f.dp(p0, q1);
}

//-----------------------------------Stormer-Verlet Method-----------------------------------------

void registerStormerVerlet(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Stormer-Verlet", [](int p){ return (TimeIntegrator*) new StormerVerletSolver(); });
}

ColVector StormerVerletSolver::oneStepSolve(TimeFunction &f, 
                                   const ColVector &x0, 
                                   const ColVector &f0, 
                                   const double &t0, 
                                   const double &step){
    auto p0 = hamiltonianGetP(x0);
    auto q0 = hamiltonianGetQ(x0);
    auto p1 = p0; auto q1 = q0;
    symEuler.oneStepP(f, p0, q0, p1, q1, step / 2);
    symEuler.oneStepQ(f, p1, q1, p0, q0, step / 2);
    return hamiltonianMerge(p0, q0);
}