#include "amgSolver.h"
#include <iostream>
using namespace std;

double u(const double &x, const double &y){
    return exp(sin(x)+y);
}

double uxx(const double &x, const double &y){
    return (-sin(x)+cos(x)*cos(x))*exp(sin(x)+y);
}

double uxy(const double &x, const double &y){
    return cos(x)*exp(sin(x)+y);
}

double uyy(const double &x, const double &y){
    return exp(sin(x)+y);
}

double f(const double &x, const double &y){
    return -(1-sin(x)+cos(x)*cos(x))*exp(sin(x)+y);
}

int IDX(const int &n, const int &i, const int &j){
    return (n-1)*(i-1)+(j-1);
}

SparseMatrix getA(const int &n){
    double h = 1.0/n;
    vector<Triple> els;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            double a, b, c;
            if(i<n/2 && j<n/2) a=1, b=0, c=1;
            if(i<n/2 && j>=n/2) a=1, b=0, c=2;
            if(i>=n/2 && j<n/2) a=2, b=0, c=1;
            if(i>=n/2 && j>=n/2) a=1, b=0.5, c=1;
            els.push_back(Triple(IDX(n,i,j), IDX(n,i,j), a*2.0/(h*h) + c*2.0/(h*h) - b/(h*h)));
            if(i>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j), -a/(h*h) + b/(2*h*h)));
            if(i<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j), -a/(h*h) + b/(2*h*h)));
            if(j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j-1), -c/(h*h) + b/(2*h*h)));
            if(j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i,j+1), -c/(h*h) + b/(2*h*h)));
            if(b && i>1 && j<n-1) els.push_back(Triple(IDX(n,i,j), IDX(n,i-1,j+1), -b/(2*h*h)));
            if(b && i<n-1 && j>1) els.push_back(Triple(IDX(n,i,j), IDX(n,i+1,j-1), -b/(2*h*h)));
        }
    return SparseMatrix((n-1)*(n-1), (n-1)*(n-1), els);
}

int main(){
    const int n = 128;
    double h = 1.0/n;
    ColVector B((n-1)*(n-1));
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            double a, b, c;
            if(i<n/2 && j<n/2) a=1, b=0, c=1;
            if(i<n/2 && j>=n/2) a=1, b=0, c=2;
            if(i>=n/2 && j<n/2) a=2, b=0, c=1;
            if(i>=n/2 && j>=n/2) a=1, b=0.5, c=1;
            B(IDX(n,i,j)) = -a*uxx(i*h, j*h) - c*uyy(i*h, j*h) + b*uxy(i*h, j*h);
            if(i==1) B(IDX(n,i,j)) += a*u(0,j*h)/(h*h) - b*u(0,j*h)/(2*h*h) + b*u(0,(j+1)*h)/(2*h*h);
            if(i==n-1) B(IDX(n,i,j)) += a*u(1,j*h)/(h*h) - b*u(1,j*h)/(2*h*h) + b*u(1,(j-1)*h)/(2*h*h);
            if(j==1) B(IDX(n,i,j)) += c*u(i*h,0)/(h*h) - b*u(i*h,0)/(2*h*h) + b*u((i+1)*h,0)/(2*h*h);
            if(j==n-1) B(IDX(n,i,j)) += c*u(i*h,1)/(h*h) - b*u(i*h,1)/(2*h*h) + b*u((i-1)*h,1)/(2*h*h);
            if(i==1 && j==n-1) B(IDX(n,i,j)) -= b*u(0,(j+1)*h)/(2*h*h);
            if(i==n-1 && j==1) B(IDX(n,i,j)) -= b*u(1,(j-1)*h)/(2*h*h);
        }
    amgSolver solver;
    solver.generateGrid(getA(n));
    ColVector sol = solver.solve(B, "FMG", 300, 1e-9);

    double maxerr = 0;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            double x = i*h, y = j*h;
            maxerr = max(maxerr, fabs(u(x,y) - sol(IDX(n,i,j))));
        }
    cout << maxerr << endl;
    return 0;
}
