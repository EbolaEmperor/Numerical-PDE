#include <bits/stdc++.h>
#include "polynomial.h"
#include "fraction.h"
using namespace std;

const int s = 3;
vector<fraction> A, b, c;

int checkB(){
    int l = 1;
    while(true){
        fraction sum(0);
        for(int j = 0; j < s; j++)
            sum += b[j] * ( c[j]^(l-1) );
        if(sum-fraction(1,l) != 0)
            return l - 1;
        l++;
    }
    return l;
}

int checkC(){
    int m = 1;
    while(true){
        for(int i = 0; i < s; i++){
            fraction sum(0);
            for(int j = 0; j < s; j++)
                sum += A[i*s+j] * ( c[j]^(m-1) );
            if(sum - (c[i]^m)*fraction(1,m) != 0)
                return m - 1;
        }
        m++;
    }
    return m;
}

int checkD(){
    int m = 1;
    while(true){
        for(int i = 0; i < s; i++){
            fraction sum(0);
            for(int j = 0; j < s; j++)
                sum += b[j] * A[j*s+i] * ( c[j]^(m-1) );
            if(sum - b[i]*fraction(1,m) + b[i]*(c[i]^m)*fraction(1,m) != 0)
                return m - 1;
        }
        m++;
    }
    return m;
}

int main(){
    c.push_back(fraction(1,4));
    c.push_back(fraction(1,2));
    c.push_back(fraction(3,4));
    vector< Polynomial<fraction> > l;
    for(int i = 0; i < 3; i++){
        Polynomial<fraction> prod = constPolynomial<fraction>(1);
        for(int j = 0; j < 3; j++){
            if(i==j) continue;
            prod = prod * (identityPolynomial<fraction>() - constPolynomial<fraction>(c[j])) / (c[i] - c[j]);
        }
        l.push_back(prod);
    }
    cout << "coefficients of A:" << endl;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            A.push_back(l[j].integral()(c[i]));
            cout << A.back() << " ";
        }
        cout << endl;
    }
    cout << "coefficients of b:" << endl;
    for(int i = 0; i < 3; i++){
        b.push_back(l[i].integral()(1));
        cout << b.back() << endl;
    }
    cout << "Properties:" << endl;
    cout << "B(" << checkB() << ")" << endl;
    cout << "C(" << checkC() << ")" << endl;
    cout << "D(" << checkD() << ")" << endl;
    return 0;
}