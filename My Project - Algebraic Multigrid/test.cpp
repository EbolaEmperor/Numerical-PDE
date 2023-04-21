#include "amgSolver.h"
#include <iostream>
using namespace std;

double sqr(double x){
    return x*x;
}

double u(double x, double y){
    return sin(x*(1-x)*y*(1-y));
}

double f(double x, double y){
    return (sqr(x-1)*sqr(x)*sqr(1-2*y) + sqr(1-2*x)*sqr(y-1)*sqr(y)) * sin((x-1)*x*(y-1)*y)
           - 2*((x-1)*x + (y-1)*y)*cos((x-1)*x*(y-1)*y);
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

SparseMatrix getA9(const int &n){
    double h = 1.0/n;
    vector<Triple> els;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            els.push_back(Triple(IDX(n,i,j), IDX(n,i,j), 8.0/(3*h*h)));
            if(i>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j), -1.0/(3*h*h)));
            if(i<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j), -1.0/(3*h*h)));
            if(j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j-1), -1.0/(3*h*h)));
            if(j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j+1), -1.0/(3*h*h)));
            if(i>1 && j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j-1), -1.0/(3*h*h)));
            if(i>1 && j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j+1), -1.0/(3*h*h)));
            if(i<n-1 && j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j-1), -1.0/(3*h*h)));
            if(i<n-1 && j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j+1), -1.0/(3*h*h)));
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
        }
    amgSolver solver;
    solver.generateGrid(getA(n));
    ColVector sol = solver.solve(b, 20, 1e-12);

    double maxerr = 0;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            double x = i*h, y = j*h;
            maxerr = max(maxerr, fabs(u(x,y) - sol(IDX(n,i,j))));
        }
    cout << maxerr << endl;
    return 0;
}