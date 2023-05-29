#include <bits/stdc++.h>
#include "IVPfactory.h"
using namespace std;

class gxTest : public TimeFunction{
private:
    double lambda, ita;
public:
    gxTest(const double &lambda, const double &ita): lambda(lambda), ita(ita) {}
    ColVector trueSol(const double &t) const{
        ColVector res(1);
        res(0) = exp(lambda*t)*(ita-1) + cos(t);
        return res;
    }
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res(1);
        res(0) = lambda * (x(0) - cos(t)) - sin(t);
        return res;
    }
};

void outputFigure(){
    gxTest g(-1e6, 1.5);
    ColVector x0(1);
    x0(0) = 1.5;
    auto solver1 = TimeIntegratorFactory.createTimeIntegrator(BDF, 1);
    solver1->solve(g, x0, 3, 30);
    auto solver2 = TimeIntegratorFactory.createTimeIntegrator(AdamsMoulton, 2);
    solver2->solve(g, x0, 3, 30);
    solver1->output("Backward Eular");
    solver2->output("Trapezoidal");
}

int main(){
    const double lambda = -1e6;
    double klist[] = {0.2, 0.1, 0.05};
    double italist[] = {1, 1.5};
    cout << "\t\tBackward Euler\tTrapezoidal" << endl;
    for(int j = 0; j < 2; j++){
        double ita = italist[j];
        cout << "eta=" << ita << "\t";
        for(int i = 0; i < 3; i++){
            if(i) cout << "\t";
            double k = klist[i];
            cout << "k=" << k << "\t";
            gxTest g(lambda, ita);
            ColVector x0(1);
            x0(0) = ita;
            auto solver1 = TimeIntegratorFactory.createTimeIntegrator(BDF, 1);
            solver1->solve(g, x0, 3, (int)round(3/k));
            cout << (g.trueSol(3) - solver1->at(3)).maxnorm() << "\t";
            auto solver2 = TimeIntegratorFactory.createTimeIntegrator(AdamsMoulton, 2);
            solver2->solve(g, x0, 3, (int)round(3/k));
            cout << (g.trueSol(3) - solver2->at(3)).maxnorm() << endl;
        }
    }
    outputFigure();
    return 0;
}
