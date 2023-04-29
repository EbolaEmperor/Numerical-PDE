#include "amgSolver.h"
#include <iostream>
#include <unordered_map>
#include <ctime>
#include <iomanip>
using namespace std;

const int n = 20;
int idx[n+1][n+1], cnt=0;

int IDX(const int &i, const int &j){
    return idx[i][j];
}

SparseMatrix genMat(){
    vector<Triple> elements;
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++){
            if(IDX(i,j)==-1) continue;
            elements.push_back(Triple(IDX(i,j), IDX(i,j), 1.0));
            if(i==n && j==n) continue;
            int m = 0;
            if(i>0 && IDX(i-1,j)!=-1) m++;
            if(i<n && IDX(i+1,j)!=-1) m++;
            if(j>0 && IDX(i,j-1)!=-1) m++;
            if(j<n && IDX(i,j+1)!=-1) m++;
            if(i>0 && IDX(i-1,j)!=-1) elements.push_back(Triple(IDX(i,j), IDX(i-1,j), -1.0/m));
            if(i<n && IDX(i+1,j)!=-1) elements.push_back(Triple(IDX(i,j), IDX(i+1,j), -1.0/m));
            if(j>0 && IDX(i,j-1)!=-1) elements.push_back(Triple(IDX(i,j), IDX(i,j-1), -1.0/m));
            if(j<n && IDX(i,j+1)!=-1) elements.push_back(Triple(IDX(i,j), IDX(i,j+1), -1.0/m));
        }
    return SparseMatrix(cnt, cnt, elements);
}

void genGraph(){
    for(int i = 0; i < n; i++){
        int u = rand()%(n+1), v = rand()%(n+1);
        if(u==0 && v==0 || u==n && v==n) continue;
        idx[u][v] = -1;
    }
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++)
            if(idx[i][j]!=-1) idx[i][j] = cnt++;
}

void GS(){
    cout << "Gauss-Seidel Method" << endl;
    int ts = clock();
    SparseMatrix A = genMat();
    ColVector x(cnt), b = ones(cnt,1);
    b(cnt-1) = 0;
    int T = 0;
    while(1){
        T++;
        ColVector nx = A.GaussSeidel(x,b);
        if(vecnorm(nx-x,0)<1e-6) break;
        x = nx;
    }
    cout << "Iteration: " << T << endl;
    cout << setprecision(3) << "Solved in " << (double)(clock()-ts)/CLOCKS_PER_SEC << "s" << endl;
    cout << setprecision(6) << "Solution: " << x(0) << endl;
    cout << "--------------------------------------------------------------" << endl;
}

void MC(){
    cout << "Mont-Carol Method" << endl;
    int ts = clock();
    int T = 1000000;
    srand(time(0));
    double step = 0;
    for(int t = 0; t < T; t++){
        int x = 0, y = 0, ans = 0;
        do{
            int m = 0;
            if(x>0 && IDX(x-1,y)!=-1) m++;
            if(x<n && IDX(x+1,y)!=-1) m++;
            if(y>0 && IDX(x,y-1)!=-1) m++;
            if(y<n && IDX(x,y+1)!=-1) m++;
            int v = rand()%m;
            ans++;
            if(x>0 && IDX(x-1,y)!=-1){
                if(v) v--;
                else { x--; continue; }
            }
            if(x<n && IDX(x+1,y)!=-1){
                if(v) v--;
                else { x++; continue; }
            }
            if(y>0 && IDX(x,y-1)!=-1){
                if(v) v--;
                else { y--; continue; }
            }
            if(y<n && IDX(x,y+1)!=-1){
                if(v) v--;
                else { y++; continue; }
            }
        }while(x!=n || y!=n);
        step += (double)ans/T;
    }
    cout << setprecision(3) << "Solved in " << (double)(clock()-ts)/CLOCKS_PER_SEC << "s" << endl;
    cout << "Solution: " << setprecision(6) << step << endl;
    cout << "--------------------------------------------------------------" << endl;
}

int main(){
    srand(19260817);
    genGraph();
    amgSolver amg;
    amg.generateGrid(genMat());
    ColVector b = ones(cnt, 1);
    b(cnt-1) = 0;
    ColVector s = amg.solve(b, "FMG", 20000, 1e-6);
    cout << "Solution: " << s(0) << endl;
    cout << "--------------------------------------------------------------" << endl;
    if(n<=50) GS();
    if(n<=5) MC();
    return 0;
}
