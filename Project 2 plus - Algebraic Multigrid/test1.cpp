#include "amgSolver.h"
#include <iostream>
using namespace std;

double u(double x, double y){
    return exp(sin(x)+y);
}

double f(double x, double y){
    return -(1-sin(x)+cos(x)*cos(x))*exp(sin(x)+y);
}

int IDX(int n, int i, int j){
    return (n-1)*(i-1)+(j-1);
}

SparseMatrix getA(const int &n){
    double h = 1.0/n;
    vector<Triple> els;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            els.push_back(Triple(IDX(n,i,j), IDX(n,i,j), 4.0/(h*h)));
            if(i>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j), -1.0/(h*h)));
            if(i<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j), -1.0/(h*h)));
            if(j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j-1), -1.0/(h*h)));
            if(j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j+1), -1.0/(h*h)));
        }
    return SparseMatrix((n-1)*(n-1), (n-1)*(n-1), els);
}

int main(){
    const int n = 256;
    double h = 1.0/n;
    ColVector b((n-1)*(n-1));
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            b(IDX(n,i,j)) = f(i*h, j*h);
            if(i==1) b(IDX(n,i,j)) += u(0, j*h) / (h*h);
            if(i==n-1) b(IDX(n,i,j)) += u(1, j*h) / (h*h);
            if(j==1) b(IDX(n,i,j)) += u(i*h, 0) / (h*h);
            if(j==n-1) b(IDX(n,i,j)) += u(i*h, 1) / (h*h);
        }
    amgSolver solver;
    solver.generateGrid(getA(n));
    ColVector sol = solver.solve(b, "FMG", 20, 1e-9);

    double maxerr = 0;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            double x = i*h, y = j*h;
            maxerr = max(maxerr, fabs(u(x,y) - sol(IDX(n,i,j))));
        }
    cout << maxerr << endl;
    return 0;
}
