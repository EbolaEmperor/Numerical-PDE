#include "matrix_fraction.h"
#include <bits/stdc++.h>
using namespace std;

void Adams_Bashforth(){
    cout << "-----------------------Adams Bashforth-----------------------" << endl;
    for(int s = 1; s <= 5; s++){
        cout << "s=" << s << "  p=" << s << "  ";
        fraction a[s+1];
        a[s] = 1;
        a[s-1] = -1;
        fracMatrix coef(s,s);
        fracMatrix rhs(s,1);
        for(int q = 1; q <= s; q++){
            for(int j = 0; j < s; j++){
                coef.element(q-1,j) = fraction(bigPow(j,q-1), bigFact(q-1));
                rhs.element(q-1,0) += fraction(bigPow(j,q), bigFact(q)) * a[j];
            }
            rhs.element(q-1,0) += fraction(bigPow(s,q), bigFact(q)) * a[s];
        }
        fracMatrix b = solve(coef, rhs);
        cout << b.reverse().T();
    }
}

void Adams_Moulton(){
    cout << "------------------------Adams Moulton------------------------" << endl;
    for(int s = 1; s <= 5; s++){
        cout << "s=" << s << "  p=" << s+1 << "  ";
        fraction a[s+1];
        a[s] = 1;
        a[s-1] = -1;
        fracMatrix coef(s+1,s+1);
        fracMatrix rhs(s+1,1);
        for(int q = 1; q <= s+1; q++){
            for(int j = 0; j <= s; j++){
                coef.element(q-1,j) = fraction(bigPow(j,q-1), bigFact(q-1));
                rhs.element(q-1,0) += fraction(bigPow(j,q), bigFact(q)) * a[j];
            }
        }
        fracMatrix b = solve(coef, rhs);
        cout << b.reverse().T();
    }
}

void BDF(){
    cout << "--------------Backward Differentiation Formula---------------" << endl;
    for(int s = 1; s <= 5; s++){
        cout << "s=" << s << "  p=" << s << "  ";
        fracMatrix coef(s+2,s+2);
        fracMatrix rhs = zeros(s+2,1);
        for(int q = 0; q <= s; q++){
            for(int j = 0; j <= s; j++){
                coef.element(q,j+1) = fraction(bigPow(j,q), bigFact(q));
            }
            if(q >= 1) coef.element(q,0) = -fraction(bigPow(s,q-1), bigFact(q-1));
        }
        coef.element(s+1,s+1) = 1;
        rhs.element(s+1,0) = 1;
        fracMatrix ab = solve(coef, rhs);
        cout << ab.reverse().T();
    }
}

int main(){
    Adams_Bashforth();
    Adams_Moulton();
    BDF();
    return 0;
}