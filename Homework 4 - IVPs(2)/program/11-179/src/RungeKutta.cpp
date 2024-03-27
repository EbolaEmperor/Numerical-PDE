#include "RungeKutta.h"
#include "NonLinearSolver.h"
#include "Polynomial.h"

//-----------------------------------Classical RK Method-----------------------------------------

void ClassicalRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    for(int i = 1; i <= n; i++){
        double t = (i-1)*timeStep;
        ColVector RK_y1 = dsol[i-1];
        ColVector RK_y2 = f(sol[i-1] + timeStep/2*RK_y1, t + timeStep/2);
        ColVector RK_y3 = f(sol[i-1] + timeStep/2*RK_y2, t + timeStep/2);
        ColVector RK_y4 = f(sol[i-1] + timeStep*RK_y3, t + timeStep);
        sol.push_back(sol[i-1] + timeStep/6 * (RK_y1 + 2*RK_y2 + 2*RK_y3 + RK_y4));
        dsol.push_back(f(sol[i],i*timeStep));
    }
}

//-------------------------------------Implicit RK Method-----------------------------------------

void ImplicitRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i],i*timeStep));
    }
}

//---------------------------------------ESDIRK Method-------------------------------------------

class ESDIRK_Function : public Function{
private:
    TimeFunction &f;
    ColVector x0;
    double b, t;
public:
    ESDIRK_Function(TimeFunction &f, const ColVector &x0, const double &b, const double &t):
        f(f), x0(x0), b(b), t(t) {}
    ColVector operator () (const ColVector &x) const{
        return f(x0 + b * x, t) - x;
    }
};

ColVector ESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    static const double aval[]={
        0,                              0,                          0,                          0,                      0,                  0,
        1.0/4.0,                        1.0/4.0,                    0,                          0,                      0,                  0,
        8611.0/62500.0,                 -1743.0/31250.0,            1.0/4.0,                    0,                      0,                  0,
        5012029.0/34652500.0,           -654441.0/2922500.0,        174375.0/388108.0,          1.0/4.0,                0,                  0,
        15267082809.0/155376265600.0,   -71443401.0/120774400.0,    730878875.0/902184768.0,    2285395.0/8070912.0,    1.0/4.0,            0,
        82889.0/524892.0,               0,                          15625.0/83664.0,            69875.0/102672.0,       -2260.0/8211.0,     1.0/4.0
    };
    static const double cval[] = {
        0,      1.0/2.0,        83.0/250.0,     31.0/50.0,      17.0/20.0,      1.0
    };
    static const Matrix A(6, 6, aval);
    static const ColVector b(6, aval+30);
    static const ColVector c(6, cval);
    static const int s = 6;
    NonLinearSolver nlsolver;

    std::vector<ColVector> y;
    y.push_back(f(U0,t0));
    for(int i = 1; i < s; i++){
        ColVector constY = U0;
        for(int j = 0; j < i; j++)
            constY = constY + step * A(i,j) * y[j];
        double t = t0 + c(i) * step;
        
        ColVector curY, nxtY = U0 + step * c(i) * y[0];
        int T = 0;
        do{
            curY = std::move(nxtY);
            nxtY = f(constY + step * A(i,i) * curY, t);
            if(++T==30) break;
        }while(relativeError(curY, nxtY) > 1e-14);
        
        if(T < 30){
            y.push_back(nxtY);
        } else {
            // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
            ESDIRK_Function esf(f, constY, step*A(i,i), t);
            y.push_back(nlsolver.solve(esf, U0 + step*c(i)*y[0]));
        }
    }

    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b(i) * y[i];
    }
    return U1;
}

//--------------------------------------Gauss-Legendre RK Methods----------------------------------------

class GaussLegendre_Function : public Function{
private:
    TimeFunction &f;
    const Matrix &A;
    const ColVector &c;
    const ColVector &U0;
    const double step, t0;
public:
    GaussLegendre_Function(TimeFunction &_f, const Matrix &_A,  const ColVector &_c, const ColVector &_U0, const double &_step, const double &_t0):
        f(_f), A(_A), c(_c), U0(_U0), step(_step), t0(_t0){}
    ColVector operator () (const ColVector &x) const{
        const int s = c.size(), m = x.size()/s;
        ColVector res(s*m);
        std::vector<ColVector> y;
        for(int i = 0; i < s; i++)
            y.push_back(x.getSubmatrix(i*m,i*m+m-1,0,0));
        for(int i = 0; i < s; i++){
            ColVector Xk = U0;
            for(int j = 0; j < s; j++)
                Xk = Xk + step * A(i,j) * y[j];
            res.setSubmatrix(i*m, 0, f(Xk, t0 + c(i)*step) - y[i]);
        }
        return res;
    }
};

GaussLegendreRKSolver::GaussLegendreRKSolver(const int &s){
    A = Matrix(s,s);
    b = ColVector(s);
    c = ColVector(s);
    
    double *coef = new double[s+1];
    for(int j = 0; j <= s; j++){
        coef[j] = sqr(fact(s))/fact(2*s) * combinational(s,j) * combinational(s+j,j);
        if((s-j)&1) coef[j] = -coef[j];
    }
    Polynomial poly(s, coef);
    auto roots = poly.roots();
    for(int i = 0; i < s; i++)
        c(i) = roots[i].real();
    c.sort();
    
    std::vector<Polynomial> l;
    for(int i = 0; i < s; i++){
        Polynomial p = constPolynomial(1);
        for(int j = 0; j < s; j++){
            if(i == j) continue;
            p = p * linearPolynomial(1, -c(j)) / (c(i) - c(j));
        }
        l.push_back(p.integral());
    }

    for(int i = 0; i < s; i++){
        b(i) = l[i](1) - l[i](0);
        for(int j = 0; j < s; j++)
            A(i,j) = l[j](c(i)) - l[j](0);
    }
}

std::vector<ColVector> GaussLegendreRKSolver_OneStepY::oneStepSolveY(TimeFunction &f, const ColVector &U0, const double &t0, const double &step, const Matrix &A, const ColVector &c){
    NonLinearSolver nlsolver;
    std::vector<ColVector> y;
    const int s = c.size(), m = U0.size();
    y.push_back(f(U0,t0));
    for(int i = 1; i < s; i++)
        y.push_back(y[0]);
    int T = 0;
    while(++T < 30){
        double maxerr = 0;
        for(int i = 0; i < s; i++){
            ColVector sumY = U0;
            for(int j = 0; j < s; j++)
                sumY = sumY + step * A(i,j) * y[j];
            ColVector nxtY = f(sumY, t0 + c(i) * step);
            maxerr = std::max(maxerr, relativeError(nxtY,y[i]));
            y[i] = nxtY;
            if(y[i].maxnorm() > 1e100){
                T = 29;
                break;
            }
        }
        if(maxerr < 1e-14) break;
    }
    if(T == 30){
        // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
        GaussLegendre_Function glf(f, A, c, U0, step, t0);
        ColVector initY(s*m), y0 = f(U0,t0);
        for(int i = 0; i < s; i++)
            initY.setSubmatrix(i*m, 0, y0);
        ColVector resY = nlsolver.solve(glf, initY);
        for(int i = 0; i < s; i++)
            y[i] = resY.getSubmatrix(i*m, i*m+m-1, 0, 0);
    }
    return y;
}

ColVector GaussLegendreRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b(i) * y[i];
    }
    return U1;
}

//----------------------------------------Embedded RK Methods----------------------------------------------

std::pair<ColVector, ColVector> EmbeddedRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    static const int s = c.size();
    std::vector<ColVector> y;
    for(int i = 0; i < s; i++){
        ColVector x1 = U0;
        for(int j = 0; j < i; j++)
            x1 = x1 + step * A(i,j) * y[j];
        y.push_back(f(x1, t0 + c(i) * step));
    }

    ColVector U1 = U0, U2 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
        U2 = U2 + step * b2(i) * y[i];
    }
    return std::make_pair(U1, U2);
}

void EmbeddedRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps){
    static const double rhomax = 5, rho = 0.8, rhomin = 0.2;
    maxTime = T;
    int gridsize = 0;
    double step = 0.1, curtime = 0.0;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    timePoint.push_back(0.0);
    while(curtime < T){
        std::pair<ColVector, ColVector> pU = oneStepSolve(f, sol[gridsize], curtime, step);
        double Eind = 0.0;
        for(int i = 0; i < pU.first.size(); i++){
            Eind += sqr( (pU.first(i) - pU.second(i)) / ( eps + fabs(sol[gridsize](i)) * eps ) );
        }
        Eind = sqrt(Eind / pU.first.size());
        if(Eind <= 1){
            sol.push_back(pU.first);
            dsol.push_back(f(pU.first, curtime));
            curtime += step;
            timePoint.push_back(curtime);
            gridsize++;
        }
        step = step * std::min( rhomax, std::max(rhomin, rho/pow(Eind,1.0/(maxorder+1))) );
    }
    std::cerr << "Grid Size: " << gridsize << std::endl;
}

//--------------------------------Fehlberg 4(5) Embedded RK Method-------------------------------------

FehlbergSolver::FehlbergSolver(){
    static const double aval[] = {
        0,              0,              0,              0,              0,              0,
        1.0/4.0,        0,              0,              0,              0,              0,
        3.0/32.0,       9.0/32.0,       0,              0,              0,              0,
        1932.0/2197.0,  -7200.0/2197.0, 7296.0/2197.0, 0,              0,              0,
        439.0/216.0,    -8.0,           3680.0/513.0,   -845.0/4104.0,  0,              0,
        -8.0/27.0,      2.0,            -3544.0/2565.0, 1859.0/4104.0,  -11.0/40.0,     0,
    };
    static const double b1val[] = {
        25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0
    };
    static const double b2val[] = {
        16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0
    };
    static const double cval[] = {
        0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0
    };
    A = Matrix(6, 6, aval);
    b1 = ColVector(6, b1val);
    b2 = ColVector(6, b2val);
    c = ColVector(6, cval);
    maxorder = 4;
}

//--------------------------------Dormand-Prince 5(4) Embedded RK Method-------------------------------------

DormandPrinceSolver::DormandPrinceSolver(){
    static const double aval[] = {
        0,              0,                  0,                0,              0,                   0,           0,
        1.0/5.0,        0,                  0,                0,              0,                   0,           0,
        3.0/40.0,       9.0/40.0,           0,                0,              0,                   0,           0,
        44.0/45.0,      -56.0/15.0,         32.0/9.0,         0,              0,                   0,           0,
        19372.0/6561.0, -25360.0/2187.0,    64448.0/6561.0,   -212.0/729.0,   0,                   0,           0,
        9017.0/3168.0,  -355.0/33.0,        46732.0/5247.0,   49.0/176.0,     -5103.0/18656.0,     0,           0,
        35.0/384.0,     0.0,                500.0/1113.0,     125.0/192.0,    -2187.0/6784.0,      11.0/84.0,   0
    };
    static const double b1val[] = {
        35.0/384.0,     0.0,                500.0/1113.0,     125.0/192.0,    -2187.0/6784.0,      11.0/84.0,   0
    };
    static const double b2val[] = {
        5179.0/57600.0,     0.0,    7571.0/16695.0,     393.0/640.0,    -92097.0/339200.0,  187.0/2100.0,   1.0/40.0
    };
    static const double cval[] = {
        0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0
    };
    A = Matrix(7, 7, aval);
    b1 = ColVector(7, b1val);
    b2 = ColVector(7, b2val);
    c = ColVector(7, cval);
    maxorder = 4;
}

//-------------------------------------Adaptive Gauss-Legendre RK Method--------------------------------------------

AdaptiveGaussLegendreRKSolver::AdaptiveGaussLegendreRKSolver(const int &s){
    GaussLegendreRKSolver gs(s);
    A = gs.getA();
    b1 = gs.getB();
    c = gs.getC();
    ColVector e(s);
    for(int i = 0; i < s-1; i++)
        e(i) = 1.0 / (i+1);
    Matrix Van(s,s);
    for(int i = 0; i < s; i++){
        double curC = 1.0;
        for(int j = 0; j < s; j++){
            Van(j,i) = curC;
            curC *= c(i);
        }
    }
    b2 = Van.solve(e);
    maxorder = 2*s;
}

std::pair<ColVector, ColVector> AdaptiveGaussLegendreRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0, U2 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
        U2 = U2 + step * b2(i) * y[i];
    }
    return std::make_pair(U1,U2);
}