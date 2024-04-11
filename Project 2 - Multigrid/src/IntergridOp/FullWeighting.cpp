#include "FullWeighting.h"

#include <cmath>

template<>
ColVector FullWeighting<1>::operator () (const ColVector &v) const{
    ColVector res(v.n/2+1);
    for(int i = 1; i < res.n-1; i++){
        res(i) = 0.25*v(2*i-1) + 0.5*v(2*i) + 0.25*v(2*i+1);
    }
    res(0) = v(0);
    res(v.n/2) = v(v.n-1);
    return res;
}

template<>
ColVector FullWeighting<2>::operator () (const ColVector &v) const{
    int n = (int)round(sqrt(v.n)) - 1;
    ColVector res( (n/2+1)*(n/2+1) );
    res(IDX(n/2,0,0)) = v(IDX(n,0,0));
    res(IDX(n/2,0,n/2)) = v(IDX(n,0,n));
    res(IDX(n/2,n/2,0)) = v(IDX(n,n,0));
    res(IDX(n/2,n/2,n/2)) = v(IDX(n,n,n));
    for(int i = 1; i < n/2; i++){
        res(IDX(n/2,i,0)) = 0.25*v(IDX(n,2*i-1,0)) + 0.5*v(IDX(n,2*i,0)) + 0.25*v(IDX(n,2*i+1,0));
        res(IDX(n/2,i,n/2)) = 0.25*v(IDX(n,2*i-1,n)) + 0.5*v(IDX(n,2*i,n)) + 0.25*v(IDX(n,2*i+1,n));
        res(IDX(n/2,0,i)) = 0.25*v(IDX(n,0,2*i-1)) + 0.5*v(IDX(n,0,2*i)) + 0.25*v(IDX(n,0,2*i+1));
        res(IDX(n/2,n/2,i)) = 0.25*v(IDX(n,n,2*i-1)) + 0.5*v(IDX(n,n,2*i)) + 0.25*v(IDX(n,n,2*i+1));
    }
    for(int i = 1; i < n/2; i++)
        for(int j = 1; j < n/2; j++){
            res(IDX(n/2,i,j)) = 0.0625 * (v(IDX(n,2*i-1,2*j-1)) + v(IDX(n,2*i-1,2*j+1)) + v(IDX(n,2*i+1,2*j-1)) + v(IDX(n,2*i+1,2*j+1)))
                                + 0.125 * (v(IDX(n,2*i-1,2*j)) + v(IDX(n,2*i+1,2*j)) + v(IDX(n,2*i,2*j-1)) + v(IDX(n,2*i,2*j+1)))
                                +  0.25 * v(IDX(n,2*i,2*j));
        }
    return res;
}

template class FullWeighting<1>;
template class FullWeighting<2>;