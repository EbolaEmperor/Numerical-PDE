#include <bits/stdc++.h>
#include "matrix_quotion.h"
#include "polynomial.h"
using namespace std;

const int s = 3;
Matrix<extendQuotion6> A(3,3), b(3,1), c(3,1);

int checkB(){
    int l = 1;
    while(true){
        extendQuotion6 sum(0);
        for(int j = 0; j < s; j++)
            sum += b(j,0) * ( c(j,0)^(l-1) );
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
            extendQuotion6 sum(0);
            for(int j = 0; j < s; j++)
                sum += A(i,j) * ( c(j,0)^(m-1) );
            if(sum - (c(i,0)^m)*fraction(1,m) != 0)
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
            extendQuotion6 sum(0);
            for(int j = 0; j < s; j++)
                sum += b(j,0) * A(j,i) * ( c(j,0)^(m-1) );
            if(sum - b(i,0)*fraction(1,m) + b(i,0)*(c(i,0)^m)*fraction(1,m) != 0)
                return m - 1;
        }
        m++;
    }
    return m;
}

void getStabilityFunction(){
    Matrix< Polynomial<extendQuotion6> > IzA(3,3), IzAeb(3,3);
    Matrix<extendQuotion6> ones(3,1), PB(1,3);
    for(int i = 0; i < s; i++){
        ones(i,0) = extendQuotion6(1);
        PB(0,i) = b(i,0);
    }
    auto eb = ones * PB;
    for(int i = 0; i < s; i++){
        IzA(i,i) = constPolynomial<extendQuotion6>(1);
        IzAeb(i,i) = constPolynomial<extendQuotion6>(1);
        for(int j = 0; j < s; j++){
            IzAeb(i,j) = IzAeb(i,j) + identityPolynomial<extendQuotion6>() * (constPolynomial<extendQuotion6>(eb(i,j)-A(i,j)));
            IzA(i,j) = IzA(i,j) + identityPolynomial<extendQuotion6>() * (constPolynomial<extendQuotion6>(-A(i,j)));
        }
    }
    cout << det(IzAeb) * constPolynomial<extendQuotion6>(60) << " / ";
    cout << det(IzA) * constPolynomial<extendQuotion6>(60) << endl;           //上下同乘60
}

int main(){
    A(0,0) = extendQuotion6( fraction(88,360), fraction(-7, 360) );
    A(0,1) = extendQuotion6( fraction(296,1800), fraction(-169, 1800) );
    A(0,2) = extendQuotion6( fraction(-2,225), fraction(3, 225) );
    A(1,0) = extendQuotion6( fraction(296,1800), fraction(169, 1800) );
    A(1,1) = extendQuotion6( fraction(88,360), fraction(7, 360) );
    A(1,2) = extendQuotion6( fraction(-2,225), fraction(-3, 225) );
    A(2,0) = extendQuotion6( fraction(16,36), fraction(-1, 36) );
    A(2,1) = extendQuotion6( fraction(16,36), fraction(1, 36) );
    A(2,2) = extendQuotion6( fraction(1,9), fraction(0, 1) );
    for(int i = 0; i < 3; i++){
        b(i,0) = A(2,i);
        for(int j = 0; j < 3; j++)
            c(i,0) += A(i,j);
    }
    cout << "B(" << checkB() << ")" << endl;
    cout << "C(" << checkC() << ")" << endl;
    cout << "D(" << checkD() << ")" << endl;
    getStabilityFunction();
    return 0;
}