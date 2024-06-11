#include "RungeKutta.h"
#include "NonLinearSolver.h"
#include "Polynomial.h"

//-----------------------------------Classical RK Method-----------------------------------------

void registerClassicalRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Classical RK", [](int p){ return (TimeIntegrator*) new ClassicalRKSolver(); });
}

ColVector ClassicalRKSolver::oneStepSolve(TimeFunction &f, const ColVector &x0, const ColVector &f0, const double &t0, const double &step){
    ColVector RK_y1 = f0;
    ColVector RK_y2 = f(x0 + step/2*RK_y1, t0 + step/2);
    ColVector RK_y3 = f(x0 + step/2*RK_y2, t0 + step/2);
    ColVector RK_y4 = f(x0 + step*RK_y3, t0 + step);
    return x0 + step/6 * (RK_y1 + 2*RK_y2 + 2*RK_y3 + RK_y4);
}

void ClassicalRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], dsol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i], i * timeStep));
    }
}

//-------------------------------------Implicit RK Method-----------------------------------------

void ImplicitRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    double errC = -(double)(1<<maxorder) / (1 - (1<<maxorder));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i],i*timeStep));
    }
}

//---------------------------------------ESDIRK Method-------------------------------------------

void registerESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("ESDIRK", [](int p){ return (TimeIntegrator*) new ESDIRKSolver(); });
}

void registerAdaptiveESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive ESDIRK", [](int p){ return (TimeIntegrator*) new AdaptiveESDIRKSolver(); });
}

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

ESDIRKSolver::ESDIRKSolver(){
    maxorder = 4;
}

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
        
        //针对线性函数的特殊优化（解热方程时有重要作用）
        if(f.isLinear()){
            ColVector _c(1);
            _c(0) = c(i);
            Matrix _a(1,1);
            _a(0,0) = A(i,i);
            y.push_back(f.solve(_a, _c, constY, t0, step));
            continue;
        }

        static const int MaxIter = 50;
        ColVector curY, nxtY = U0 + step * c(i) * y[0];
        int T = 0;
        do{
            curY = std::move(nxtY);
            nxtY = f(constY + step * A(i,i) * curY, t);
            if(++T==MaxIter) break;
            if(nxtY.maxnorm() > 1e100){
                T = MaxIter;
                break;
            }
        }while(relativeError(curY, nxtY) > 1e-14);
        
        if(T < MaxIter){
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

AdaptiveESDIRKSolver::AdaptiveESDIRKSolver(){
    maxorder = 4;
}

std::pair<ColVector, ColVector> AdaptiveESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    ColVector U1 = esdirk.oneStepSolve(f, U0, t0, step);
    ColVector U2_mid = esdirk.oneStepSolve(f, U0, t0, step/2);
    ColVector U2 = esdirk.oneStepSolve(f, U2_mid, t0+step/2, step/2);
    return std::make_pair(U1,U2);
}

//--------------------------------------Embedded ESDIRK Methods----------------------------------------

void registerEmbeddedESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Embedded ESDIRK", [](int p){ return (TimeIntegrator*) new EmbeddedESDIRKSolver(p); });
}

EmbeddedESDIRKSolver::EmbeddedESDIRKSolver(const int &p){
    if(p < 4 || p > 5){
        std::cerr << "The order of an Embedded ESDIRK method could only be 4 or 5." << std::endl;
        exit(-1);
    }
    if(p == 4){
        static const double aval[]={
            0,                              0,                          0,                          0,                      0,                  0,                0,
            1.0/8.0,                     1.0/8.0,                 0,                          0,                      0,                  0,                 0,
            -39188347878.0/1513744654945.0,  -39188347878.0/1513744654945.0,  1.0/8.0,         0,                      0,                  0,                  0,
            1748874742213.0/5168247530883.0,  1748874742213.0/5168247530883.0,  -1748874742213.0/5795261096931.0,          1.0/8.0,                0,                  0,                 0,
            -6429340993097.0/17896796106705.0,  -6429340993097.0/17896796106705.0,  9711656375562.0/10370074603625.0,    1137589605079.0/3216875020685.0,    1.0/8.0,            0,                   0,
            405169606099.0/1734380148729.0,  405169606099.0/1734380148729.0,  -264468840649.0/6105657584947.0,     118647369377.0/6233854714037.0,       683008737625.0/4934655825458.0,     1.0/8.0,              0,
            -5649241495537.0/14093099002237.0,  -5649241495537.0/14093099002237.0,  5718691255176.0/6089204655961.0,  2199600963556.0/4241893152925.0,  8860614275765.0/11425531467341.0,  -3696041814078.0/6641566663007.0,       1.0/8.0
        };
        static const double b2val[] = {
            -1517409284625.0/6267517876163.0,  -1517409284625.0/6267517876163.0,      8291371032348.0/12587291883523.0,        5328310281212.0/10646448185159.0,     5405006853541.0/7104492075037.0,      -4254786582061.0/7445269677723.0,      19.0/140.0
        };
        static const double cval[] = {
            0,      1.0/4.0,        (2.0-sqrt(2.0))/8.0,     1.0/2.0,      395.0/567.0,      89.0/126.0,       1.0
        };
        A = Matrix(7, 7, aval);
        b1 = ColVector(7, aval+42);
        b2 = ColVector(7, b2val);
        c = ColVector(7, cval);
        maxorder = 4;
    } else {
        static const double aval[]={
            0,                              0,                          0,                       0,                                                 0,                          0,                                    0,                              0,
            1.0/7.0,                     1.0/7.0,                       0,                       0,                                                 0,                          0,                                    0,                              0,
            1521428834970.0/8822750406821.0,  1521428834970.0/8822750406821.0,  1.0/7.0,         0,                                                 0,                          0,                                    0,                              0,
            5338711108027.0/29869763600956.0,  5338711108027.0/29869763600956.0,  1483184435021.0/6216373359362.0,          1.0/7.0,                0,                          0,                                    0,                              0,
            2264935805846.0/12599242299355.0,  2264935805846.0/12599242299355.0,  1330937762090.0/13140498839569.0,    -287786842865.0/17211061626069.0,    1.0/7.0,            0,                                    0,                              0,
            118352937080.0/527276862197.0,  118352937080.0/527276862197.0,  -2960446233093.0/7419588050389.0,     -3064256220847.0/46575910191280.0,       6010467311487.0/7886573591137.0,     1.0/7.0,              0,                              0,
            1134270183919.0/9703695183946.0,  1134270183919.0/9703695183946.0,  4862384331311.0/10104465681802.0,  1127469817207.0/2459314315538.0,  -9518066423555.0/11243131997224.0,  -811155580665.0/7490894181109.0,       1.0/7.0,              0,
            2162042939093.0/22873479087181.0,  2162042939093.0/22873479087181.0,  -4222515349147.0/9397994281350.0,  3431955516634.0/4748630552535.0,  -374165068070.0/9085231819471.0,  -1847934966618.0/8254951855109.0,  5186241678079.0/7861334770480.0,       1.0/7.0
        };
        static const double b2val[] = {
            701879993119.0/7084679725724.0,  701879993119.0/7084679725724.0,  -8461269287478.0/14654112271769.0,      6612459227430.0/11388259134383.0,        2632441606103.0/12598871370240.0,     -2147694411931.0/10286892713802.0,      4103061625716.0/6371697724583.0,      36.0/233.0
        };
        static const double cval[] = {
            0,      2.0/7.0,        (2.0+sqrt(2.0))/7.0,     150.0/203.0,      27.0/46.0,      473.0/532.0,       30.0/83.0,       1.0
        };
        A = Matrix(8, 8, aval);
        b1 = ColVector(8, aval+56);
        b2 = ColVector(8, b2val);
        c = ColVector(8, cval);
        maxorder = 5;
    }
}

std::pair<ColVector, ColVector> EmbeddedESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    NonLinearSolver nlsolver;
    std::vector<ColVector> y;
    const int s = c.size();
    y.push_back(f(U0,t0));

    for(int i = 1; i < s; i++){
        ColVector constY = U0;
        for(int j = 0; j < i; j++)
            constY = constY + step * A(i,j) * y[j];
        double t = t0 + c(i) * step;
        
        //针对线性函数的特殊优化（解热方程时有重要作用）
        if(f.isLinear()){
            ColVector _c(1);
            _c(0) = c(i);
            Matrix _a(1,1);
            _a(0,0) = A(i,i);
            y.push_back(f.solve(_a, _c, constY, t0, step));
            continue;
        }

        static const int MaxIter = 50;
        ColVector curY, nxtY = U0 + step * c(i) * y[0];
        int T = 0;
        do{
            curY = std::move(nxtY);
            nxtY = f(constY + step * A(i,i) * curY, t);
            if(++T==MaxIter) break;
            if(nxtY.maxnorm() > 1e100){
                T = MaxIter;
                break;
            }
        }while(relativeError(curY, nxtY) > 1e-14);
        
        if(T < MaxIter){
            y.push_back(nxtY);
        } else {
            // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
            ESDIRK_Function esf(f, constY, step*A(i,i), t);
            y.push_back(nlsolver.solve(esf, U0 + step*c(i)*y[0]));
        }
    }

    ColVector U1 = U0, U2 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
        U2 = U2 + step * b2(i) * y[i];
    }
    return std::make_pair(U1,U2);
}

//--------------------------------------Gauss-Legendre RK Methods----------------------------------------

void registerGaussLegendre(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Gauss-Legendre", [](int p){ return (TimeIntegrator*) new GaussLegendreRKSolver(p); });
}

GaussLegendreRKSolver::GaussLegendreRKSolver(const int &s){
    if(s <= 0){
        std::cerr << "[Error] The stage of a Gauss-Legendre solver should be at least 1." << std::endl;
        exit(-1);
    }
    A = Matrix(s,s);
    b = ColVector(s);
    c = ColVector(s);
    
    double *coef = new double[s+1];
    for(int j = 0; j <= s; j++){
        coef[j] = sqr(fact(s))/fact(2*s) * binom(s,j) * binom(s+j,j);
        if((s-j)&1) coef[j] = -coef[j];
    }
    Polynomial poly(s, coef);
    auto roots = poly.roots();
    for(int i = 0; i < s; i++)
        c(i) = roots[i].real();
    c.sort();
    computeAandB();
    maxorder = 2*s;
}

//--------------------------------------Radau IIA RK Methods----------------------------------------

void registerRadauIIA(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Radau-IIA", [](int p){ return (TimeIntegrator*) new RadauIIARKSolver(p); });
}

RadauIIARKSolver::RadauIIARKSolver(const int &s){
    if(s <= 0){
        std::cerr << "[Error] The stage of a Radau-IIA solver should be at least 1." << std::endl;
        exit(-1);
    }
    A = Matrix(s,s);
    b = ColVector(s);
    c = ColVector(s);

    Polynomial poly = pow(linearPolynomial(1,0),s-1) * pow(linearPolynomial(1,-1),s);
    for(int i = 0; i < s-1; i++) poly = poly.derivative();
    auto p = poly.roots();
    for(int i = 0; i < p.size(); i++)
        c(i) = p[i].real();
    c.sort();
    computeAandB();
    maxorder = 2*s-1;
}

//-------------------------Some Generic Functions for Adaptive RK Methods------------------------------

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

void AdaptiveRKSolver::solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps){
    static const double rhomax = 5, rhomin = 0.2, rho = 0.2;
    maxTime = T;
    if(maxStep==-1) maxStep = T/100;
    int gridsize = 0;
    double step = T/1e4, curtime = 0.0;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    timePoint.push_back(0.0);
    while(curtime < T){
        if(curtime + step >= T) step = T - curtime;
        std::pair<ColVector, ColVector> pU = oneStepSolve(f, sol[gridsize], curtime, step);
        double Eind = 0.0;
        for(int i = 0; i < pU.first.size(); i++){
            Eind += sqr( (pU.first(i) - pU.second(i)) / ( eps + fabs(sol[gridsize](i)) * eps ) );
        }
        Eind = sqrt(Eind / pU.first.size());
        if(Eind <= 1){
            curtime += step;
            sol.push_back(pU.first);
            dsol.push_back(f(pU.first, curtime));
            timePoint.push_back(curtime);
            gridsize++;
        }
        step = step * std::min( rhomax, std::max(rhomin, pow(rho/Eind, 1.0/(maxorder+1))) );
        if(step > maxStep) step = maxStep;
    }
}

//--------------------------------Fehlberg 4(5) Embedded RK Method-------------------------------------

void registerFehlberg(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Fehlberg", [](int p){ return (TimeIntegrator*) new FehlbergSolver(); });
}

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

//--------------------------------Dormand-Prince 8(7) Embedded RK Method-------------------------------------

void registerDormandPrince8(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Dormand-Prince 8(7)", [](int p){ return (TimeIntegrator*) new DormandPrince8Solver(); });
}

DormandPrince8Solver::DormandPrince8Solver(){
    static const double aval[] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0/18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0/48.0, 1.0/16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0/32.0, 0.0, 3.0/32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0, -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0, 800635310.0/3783071287.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0, -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0, 393006217.0/1396673457.0, 123872331.0/1001029789.0, 0.0, 0.0, 0.0, 0.0,
        -1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0, 1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0, -45442868181.0/3398467696.0, 3065993473.0/597172653.0, 0.0, 0.0, 0.0,
        185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0, -477755414.0/1098053517.0, -703635378.0/230739211.0, 5731566787.0/1027545527.0, 5232866602.0/850066563.0, -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0, 0.0,0.0,
        403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0, -411421997.0/543043805.0, 652783627.0/914296604.0, 11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0, 248638103.0/1413531060.0, 0.0, 0.0,
    };
    static const double b1val[] = {
        14005451.0/335480064.0, 0.0, 0.0, 0.0, 0.0, -59238493.0/1068277825.0, 181606767.0/758867731.0, 561292985.0/797845732.0, -1041891430.0/1371343529.0, 760417239.0/1151165299.0, 118820643.0/751138087.0, -528747749.0/2220607170.0, 1.0/4.0
    };
    static const double b2val[] = {
        13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0, -808719846.0/976000145.0, 1757004468.0/5645159321.0, 656045339.0/265891186.0, -3867574721.0/1518517206.0, 465885868.0/322736535.0, 53011238.0/667516719.0, 2.0/45.0, 0.0
    };
    static const double cval[] = {
        0.0, 1.0/18.0, 1.0/12.0, 1.0/8.0, 5.0/16.0, 3.0/8.0, 59.0/400.0, 93.0/200.0, 5490023248.0/9719169821.0, 13.0/20.0, 1201146811.0/1299019798.0, 1.0, 1.0
    };
    A = Matrix(13, 13, aval);
    b1 = ColVector(13, b1val);
    b2 = ColVector(13, b2val);
    c = ColVector(13, cval);
    maxorder = 7;
}

//--------------------------------Dormand-Prince 5(4) Embedded RK Method-------------------------------------

void registerDormandPrince(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Dormand-Prince", [](int p){ return (TimeIntegrator*) new DormandPrinceSolver(); });
}

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

void registerAdaptiveGaussLegendre(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive Gauss-Legendre", [](int p){ return (TimeIntegrator*) new AdaptiveGaussLegendreRKSolver(p); });
}

AdaptiveGaussLegendreRKSolver::AdaptiveGaussLegendreRKSolver(const int &s){
    maxorder = 2*s;
    CollocationRKSolver *p = new GaussLegendreRKSolver(s);
    computeCoef(p);
}

//-------------------------------------Adaptive Radau IIA RK Method--------------------------------------------

void registerAdaptiveRadauIIA(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive Radau-IIA", [](int p){ return (TimeIntegrator*) new AdaptiveRadauIIARKSolver(p); });
}

AdaptiveRadauIIARKSolver::AdaptiveRadauIIARKSolver(const int &s){
    maxorder = 2*s - 1;
    CollocationRKSolver *p = new RadauIIARKSolver(s);
    computeCoef(p);
}

//-------------------------Some Generic Functions for Collocation Methods---------------------------

void CollocationRKSolver::computeAandB(){
    const int s = c.size();
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

class CollocationSolver_Function : public Function{
private:
    TimeFunction &f;
    const Matrix &A;
    const ColVector &c;
    const ColVector &U0;
    const double step, t0;
public:
    CollocationSolver_Function(TimeFunction &_f, const Matrix &_A,  const ColVector &_c, const ColVector &_U0, const double &_step, const double &_t0):
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

std::vector<ColVector> CollocationRKSolver_OneStepY::oneStepSolveY(TimeFunction &f, const ColVector &U0, const double &t0, const double &step, const Matrix &A, const ColVector &c){
    NonLinearSolver nlsolver;
    std::vector<ColVector> y;
    const int s = c.size(), m = U0.size();
    y.push_back(f(U0,t0));
    for(int i = 1; i < s; i++)
        y.push_back(y[0]);
    int T = 0;
    static const int MaxIter = 50;
    while(!f.isLinear() && ++T < MaxIter){
        double maxerr = 0;
        for(int i = 0; i < s; i++){
            ColVector sumY = U0;
            for(int j = 0; j < s; j++)
                sumY = sumY + step * A(i,j) * y[j];
            ColVector nxtY = f(sumY, t0 + c(i) * step);
            maxerr = std::max(maxerr, relativeError(nxtY,y[i]));
            y[i] = nxtY;
            if(y[i].maxnorm() > 1e100){
                T = MaxIter-1;
                break;
            }
        }
        if(maxerr < 1e-14) break;
    }
    if(f.isLinear() || T == MaxIter){
        // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
        ColVector resY;
        if(!f.isLinear()){
            CollocationSolver_Function glf(f, A, c, U0, step, t0);
            ColVector initY(s*m), y0 = f(U0,t0);
            for(int i = 0; i < s; i++)
                initY.setSubmatrix(i*m, 0, y0);
            resY = nlsolver.solve(glf, initY);
        } else {
            // 对于线性函数f，可以直接求解线性方程组，不需要迭代，这将大大优化求解速度
            resY = f.solve(A, c, U0, t0, step);
        }
        for(int i = 0; i < s; i++)
            y[i] = resY.getSubmatrix(i*m, i*m+m-1, 0, 0);
    }
    return y;
}

ColVector ConstStepCollocationRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b(i) * y[i];
    }
    return U1;
}

ColVector AdaptiveCollocationRKSolver::trueOneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
    }
    return U1;
}

std::pair<ColVector, ColVector> AdaptiveCollocationRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const double &t0, const double &step){
    ColVector U1 = trueOneStepSolve(f, U0, t0, step);
    ColVector U2_mid = trueOneStepSolve(f, U0, t0, step/2);
    ColVector U2 = trueOneStepSolve(f, U2_mid, t0+step/2, step/2);
    return std::make_pair(U1,U2);
}

void AdaptiveCollocationRKSolver::computeCoef(CollocationRKSolver *p){
    A = p->getA();
    b1 = p->getB();
    c = p->getC();
}
