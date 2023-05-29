#include <bits/stdc++.h>
using namespace std;

const double pi = acos(-1);

double initial(const double &x){
    if(0.45 <= x && x < 0.5) return 20*(x-0.45);
    else if(0.5 <= x && x < 0.55) return -20*(x-0.55);
    else return 0;
}

double trueSol(const double &x, const double &t){
    double sum = 0;
    for(int k = 1; k <= 20; k++){
        double Ak = 40/(k*k*pi*pi) * (-sin(9*k*pi/20) + 2*sin(k*pi/2) - sin(11*k*pi/20));
        sum += Ak * exp(-k*k*pi*pi*t) * sin(k*pi*x);
    }
    return sum;
}

int main(int argc, char* argv[]){
    double t = stod(argv[1]);
    ofstream fout("result.txt");
    for(int i = 0; i <= 1000; i++){
        double val = t ? trueSol(0.001*i, t) : initial(0.001*i);
        fout << 0.001*i << " " << val << "\n";
    }
    fout.close();
    return 0;
}