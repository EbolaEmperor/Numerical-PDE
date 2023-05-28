#include "LMMs.h"
#include "NonLinearSolver.h"
#include "RungeKutta.h"

//-------------------------------------LMM Solver---------------------------------------

void LMMSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(sol[0],0));
    if(s >= 2){
        TimeIntegrator* initSolver = new RadauIIARKSolver(p/2+1);
        initSolver->solve( f, x0, (s-1)*timeStep, s-1);
        for(int i = 1; i < s; i++){
            sol.push_back(initSolver->at(i*timeStep));
            dsol.push_back(f(sol[i], i*timeStep));
        }
        delete initSolver;
    }
    for(int i = s; i <= n; i++){
        sol.push_back(oneStepSolve(f, i));
        dsol.push_back(f(sol[i], i*timeStep));
    }
}

//-----------------------------------Adams Bashforth--------------------------------------

void registerAdamsBashforth(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adams-Bashforth", [](int p){ return (TimeIntegrator*) new AdamsBashforthSolver(p); });
}

AdamsBashforthSolver::AdamsBashforthSolver(const int &_p){
    if(_p <= 0){
        std::cerr << "[Error] The order of an Adams-Bashforth solver should be at least 1." << std::endl;
        exit(-1);
    }
    method = "Adams-Bashforth";
    s = p = _p;
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

ColVector AdamsBashforthSolver::oneStepSolve(TimeFunction &f, const int &i){
    ColVector nxtsol = sol[i-1];
    for(int j = 0; j < s; j++)
        nxtsol = nxtsol + timeStep * beta(j) * dsol[i-s+j];
    return nxtsol;
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

void registerAdamsMoulton(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adams-Moulton", [](int p){ return (TimeIntegrator*) new AdamsMoultonSolver(p); });
}

AdamsMoultonSolver::AdamsMoultonSolver(const int &_p){
    if(_p <= 1){
        std::cerr << "[Error] The order of an Adams-Moulton solver should be at least 2." << std::endl;
        exit(-1);
    }
    method = "Adams-Moulton";
    p = _p;
    s = _p - 1;
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

ColVector AdamsMoultonSolver::oneStepSolve(TimeFunction &f, const int &i){
    ColVector C = sol[i-1];
    for(int j = 0; j < s; j++)
        C = C + timeStep * beta(j) * dsol[i-s+j];
    ColVector nxtsol = C + timeStep * beta(s) * f(sol[i-1], i*timeStep), cursol;
    int times = 0;
    do{
        times++;
        cursol = std::move(nxtsol);
        nxtsol = C + timeStep * beta(s) * f(cursol, i*timeStep);
        if(nxtsol.maxnorm() > 1e100){
            times = 40;
            break;
        }
    } while( relativeError(cursol, nxtsol) > 1e-15 && times < 40);
    if(times == 40){
        //当不动点迭代不收敛时，调用非线性求解器
        static NonLinearSolver nlsolver;
        Function_ImplicitLMM flmm(f, timeStep * beta(s), i * timeStep, C);
        nxtsol = nlsolver.solve(flmm, sol[i-1]);
    }
    return nxtsol;
}

//------------------------Backward Differential Formula-----------------------------

void registerBDF(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("BDF", [](int p){ return (TimeIntegrator*) new BDFSolver(p); });
}

BDFSolver::BDFSolver(const int &_p){
    if(_p <= 0){
        std::cerr << "[Error] The order of an BDF solver should be at least 1." << std::endl;
        exit(-1);
    }
    method = "BDF";
    p = s = _p;
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

ColVector BDFSolver::oneStepSolve(TimeFunction &f, const int &i){
    ColVector C(sol[0].size());
    for(int j = 0; j < s; j++)
        C = C - alpha(j) * sol[i-s+j];
    ColVector nxtsol = C + timeStep * beta * f(sol[i-1], i*timeStep), cursol;
    int times = 0;
    do{
        times++;
        cursol = std::move(nxtsol);
        nxtsol = C + timeStep * beta * f(cursol, i*timeStep);
        if(nxtsol.maxnorm() > 1e100){
            times = 40;
            break;
        }
    } while(relativeError(cursol, nxtsol) > 1e-15 && times < 40);
    if(times == 40){
        //当不动点迭代不收敛时，调用非线性求解器
        static NonLinearSolver nlsolver;
        Function_ImplicitLMM flmm(f, timeStep * beta, i * timeStep, C);
        nxtsol = nlsolver.solve(flmm, sol[i-1]);
    }
    return nxtsol;
}