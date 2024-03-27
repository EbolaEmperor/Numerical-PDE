#include "LMMs.h"
#include "NonLinearSolver.h"

//-----------------------------------Adams Bashforth--------------------------------------

AdamsBashforthSolver::AdamsBashforthSolver(const int &p){
    Matrix coef(p,p);
    ColVector rhs(p);
    for(int q = 1; q <= p; q++){
        for(int j = 0; j < p; j++){
            coef(q-1,j) = pow(j,q-1)/fact(q-1);
        }
        rhs(q-1) = pow(p,q)/fact(q) - pow(p-1,q)/fact(q);
    }
    beta = coef.solve(rhs);
}

void AdamsBashforthSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    int s = beta.size();
    if(s >= 2){
        TimeIntegrator* initSolver = new AdamsBashforthSolver(1);
        initSolver->solve(f, x0, (s-1)*timeStep, s);
        for(int i = 1; i < s; i++)
           sol.push_back(initSolver->at(i*timeStep));
    }
    for(int i = 0; i < s; i++)
        dsol.push_back(f(sol[i], i*timeStep));
    for(int i = s; i <= n; i++){
        sol.push_back(sol[i-1]);
        for(int j = 0; j < s; j++)
            sol[i] = sol[i] + timeStep * beta(j) * dsol[i-s+j];
        dsol.push_back(f(sol[i], i*timeStep));
    }
    
}


//-----------------------------Non-Linear Solver Function for Implicit LMM--------------------------------

class Function_ImplicitLMM : public Function{
private:
    const double a, t;
    const ColVector &C;
    TimeFunction &f;
public:
    Function_ImplicitLMM(TimeFunction &_f, const double &_a, const double &_t, const ColVector &_C):
        f(_f), a(_a), t(_t), C(_C) {}
    ColVector operator () (const ColVector &x) const{
        return C + a * f(x, t) - x;
    }
};


//--------------------------------------Adams Moulton------------------------------------------

AdamsMoultonSolver::AdamsMoultonSolver(const int &p){
    if(p <= 1){
        std::cerr << "[Error] The order of an Adams-Moulton solver should be at least 2." << std::endl;
        exit(-1);
    }
    Matrix coef(p,p);
    ColVector rhs(p);
    for(int q = 1; q <= p; q++){
        for(int j = 0; j < p; j++){
            coef(q-1,j) = pow(j,q-1)/fact(q-1);
        }
        rhs(q-1) = pow(p-1,q)/fact(q) - pow(p-2,q)/fact(q);
    }
    beta = coef.solve(rhs);
}

void AdamsMoultonSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    int s = beta.size() - 1;
    if(s >= 2){
        TimeIntegrator* initSolver = new AdamsMoultonSolver(2);
        initSolver->solve(f, x0, (s-1)*timeStep, s);
        for(int i = 1; i < s; i++)
           sol.push_back(initSolver->at(i*timeStep));
    }
    for(int i = 0; i < s; i++)
        dsol.push_back(f(sol[i], i*timeStep));
    for(int i = s; i <= n; i++){
        ColVector C = sol[i-1];
        for(int j = 0; j < s; j++)
            C = C + timeStep * beta(j) * dsol[i-s+j];
        ColVector nxtsol = C + timeStep * beta(s) * f(sol[i-1], i*timeStep), cursol;
        int times = 0;
        do{
            times++;
            cursol = std::move(nxtsol);
            nxtsol = C + timeStep * beta(s) * f(cursol, i*timeStep);
        } while(relativeError(cursol, nxtsol) > 1e-15 && times < 40);
        if(times == 40){
            //当不动点迭代不收敛时，调用非线性求解器
            static NonLinearSolver nlsolver;
            Function_ImplicitLMM flmm(f, timeStep * beta(s), i * timeStep, C);
            nxtsol = nlsolver.solve(flmm, sol[i-1]);
        }
        sol.push_back(nxtsol);
        dsol.push_back(f(nxtsol, i*timeStep));
    }
}

//------------------------Backward Differential Formula-----------------------------

BDFSolver::BDFSolver(const int &p){
    Matrix coef(p+2,p+2);
    ColVector rhs(p+2);
    for(int q = 0; q <= p; q++){
        for(int j = 0; j <= p; j++){
            coef(q,j+1) = pow(j,q)/fact(q);
        }
        if(q >= 1)
            coef(q,0) = -pow(p,q-1)/fact(q-1);
    }
    coef(p+1,p+1) = 1;
    rhs(p+1) = 1;
    ColVector ab = coef.solve(rhs);
    beta = ab(0);
    alpha = ab.getSubmatrix(1,p,0,0);
}

void BDFSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    int s = alpha.size();
    if(s >= 2){
        TimeIntegrator* initSolver = new BDFSolver(1);
        initSolver->solve(f, x0, (s-1)*timeStep, s);
        for(int i = 1; i < s; i++)
           sol.push_back(initSolver->at(i*timeStep));
    }
    for(int i = 0; i < s; i++)
        dsol.push_back(f(sol[i], i*timeStep));
    for(int i = s; i <= n; i++){
        ColVector C(x0.size());
        for(int j = 0; j < s; j++)
            C = C - alpha(j) * sol[i-s+j];
        ColVector nxtsol = C + timeStep * beta * f(sol[i-1], i*timeStep), cursol;
        int times = 0;
        do{
            times++;
            cursol = std::move(nxtsol);
            nxtsol = C + timeStep * beta * f(cursol, i*timeStep);
        } while(relativeError(cursol, nxtsol) > 1e-15 && times < 40);
        if(times == 40){
            //当不动点迭代不收敛时，调用非线性求解器
            static NonLinearSolver nlsolver;
            Function_ImplicitLMM flmm(f, timeStep * beta, i * timeStep, C);
            nxtsol = nlsolver.solve(flmm, sol[i-1]);
        }
        sol.push_back(nxtsol);
        dsol.push_back(f(nxtsol, i*timeStep));
    }
}