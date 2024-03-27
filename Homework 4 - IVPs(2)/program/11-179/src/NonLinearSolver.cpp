#include "NonLinearSolver.h"
#include <cmath>
using namespace std;

//----------------------------------Function-----------------------------------

Matrix Function::jacobi(const ColVector &x) const{
    ColVector y0 = (*this)(x);
    const int n = x.size();
    const int m = y0.size();
    Matrix h(m,n);
    for(int j = 0; j < n; j++){
        ColVector d(n);
        double scale = max(1.0, min((*this)(x).maxnorm(), 1e7));
        d(j) = 1e-13 * scale;
        double alpha = 5e12 / scale;
        ColVector g1 = alpha*( (*this)(x+d) - (*this)(x-d) );
        h.setSubmatrix(0,j,g1);
    }
    return h;
}

double Function::fval(const ColVector &x) const{
    ColVector y = (*this)(x);
    return 0.5*(y.T()*y);
}

ColVector Function::gradval(const ColVector &x) const{
    return jacobi(x).T() * (*this)(x);
}

//------------------------Non-Linear Solver (Newton Method)------------------------------------

ColVector NonLinearSolver::solve(Function &F, ColVector current){
    return solve(F, current, 1e-15);
}

ColVector NonLinearSolver::solve(Function &F, ColVector current, const double &err){
    return solve(F, current, err, 100);
}

ColVector NonLinearSolver::solve(Function &F, ColVector current, const double &err, const int &MAXN){
    const int n = current.size();
    Matrix J = F.jacobi(current);
    ColVector v = F(current);
    int T = 0;
    while(v.vecnorm(2) >= err && ++T <= MAXN){
        if(det(J) == 0) J = J + 0.01 * eye(v.size());
        ColVector d = J.solve(v);
        current = current - d;
        J = F.jacobi(current);
        v = F(current);
    }
    return current;
/************************* (Non-Linear FR Method) ****************************
    const int n = current.size();
    ColVector g = F.gradval(current), d = -g;
    int step = 0;
    while(g.vecnorm(2)>=err && F.fval(current)>=err && ++step <= MAXN){
        double alpha = search(F, current, d, g);
        current = current + alpha * d;
        if(step % (n+1) == 0){
            g = F.gradval(current);
            d = -g;
        } else {
            double tmp = g.sqrsum();
            g = F.gradval(current);
            double beta = g.sqrsum() / tmp;
            d = -g + beta * d;
        }
    }
    return current;
*******************************************************************************/
}

double NonLinearSolver::search(Function &f, const ColVector &initial, const ColVector &direct, const ColVector &grad){
    int step = 0;
    static const double sigma = 0.4, rho = 0.55;
    double alpha = 1, f0 = f.fval(initial), g0 = direct.T()*grad;
    while(++step < 10){
        double fa = f.fval(initial + alpha * direct);
        if(fa <= f0 + sigma * alpha * g0)
            return alpha;
        alpha = alpha * rho;
    }
    return step==10 ? 1 : alpha;
}
