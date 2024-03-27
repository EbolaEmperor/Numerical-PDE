#include <bits/stdc++.h>
#include "fraction.h"
using namespace std;

const int s = 4;
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

int check4Order(){
    vector<fraction> sum(4);
    sum[0] = 0;
    for(int j = 0; j < s; j++)
        sum[0] += b[j];
    if(sum[0] != 1) return 0;

    sum[0] = 0;
    for(int j = 0; j < s; j++)
        for(int k = 0; k < s; k++)
            sum[0] += b[j] * A[j*s+k];
    if(sum[0]*2 != 1) return 1;

    sum[0] = sum[1] = 0;
    for(int j = 0; j < s; j++)
        for(int k = 0; k < s; k++)
            for(int l = 0; l < s; l++){
                sum[0] += b[j] * A[j*s+k] * A[j*s+l];
                sum[1] += b[j] * A[j*s+k] * A[k*s+l];
            }
    if(sum[0]*3 != 1 || sum[1]*6 != 1) return 2;

    sum[0] = sum[1] = sum[2] = sum[3] = 0;
    for(int j = 0; j < s; j++)
        for(int k = 0; k < s; k++)
            for(int l = 0; l < s; l++)
                for(int m = 0; m < s; m++){
                    sum[0] += b[j] * A[j*s+k] * A[j*s+l] * A[j*s+m];
                    sum[1] += b[j] * A[j*s+k] * A[k*s+l] * A[j*s+m];
                    sum[2] += b[j] * A[j*s+k] * A[k*s+l] * A[k*s+m];
                    sum[3] += b[j] * A[j*s+k] * A[k*s+l] * A[l*s+m];
                }
    if(sum[0]*4 != 1 || sum[1]*8 != 1 || sum[2]*12 != 1 || sum[3]*24 != 1) return 3;
    return 4;
}

int main(){
    b.push_back(fraction(1,6));
    b.push_back(fraction(1,3));
    b.push_back(fraction(1,3));
    b.push_back(fraction(1,6));
    c.push_back(fraction(0,1));
    c.push_back(fraction(1,2));
    c.push_back(fraction(1,2));
    c.push_back(fraction(1,1));
    A.resize(16);
    A[4] = fraction(1,2);
    A[9] = fraction(1,2);
    A[14] = fraction(1,1);
    cout << "Order: >= " << check4Order() << endl;
    cout << "B(" << checkB() << ")" << endl;
    return 0;
}