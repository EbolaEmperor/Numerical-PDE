#include "LinearInterpolation.h"

#include <cmath>

template<>
ColVector LinearInterpolation<1>::operator () (const ColVector &v) const{
    ColVector res(v.n*2-1);
    for(int i = 0; i < v.n; i++){
        res(2*i) = v(i);
        if(i) res(2*i-1) += 0.5*v(i);
        if(i<v.n-1) res(2*i+1) += 0.5*v(i);
    }
    return res;
}

template<>
ColVector LinearInterpolation<2>::operator () (const ColVector &v) const{
    int n = (int)round(sqrt(v.n)) - 1;
    ColVector res( (2*n+1)*(2*n+1) );
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++){
            res(IDX(2*n,2*i,2*j)) = v(IDX(n,i,j));
            if(i) res(IDX(2*n,2*i-1,2*j)) += 0.5*v(IDX(n,i,j));
            if(i<n) res(IDX(2*n,2*i+1,2*j)) += 0.5*v(IDX(n,i,j));
            if(j) res(IDX(2*n,2*i,2*j-1)) += 0.5*v(IDX(n,i,j));
            if(j<n) res(IDX(2*n,2*i,2*j+1)) += 0.5*v(IDX(n,i,j));
            if(i && j) res(IDX(2*n,2*i-1,2*j-1)) += 0.25*v(IDX(n,i,j));
            if(i && j<n) res(IDX(2*n,2*i-1,2*j+1)) += 0.25*v(IDX(n,i,j));
            if(i<n && j) res(IDX(2*n,2*i+1,2*j-1)) += 0.25*v(IDX(n,i,j));
            if(i<n && j<n) res(IDX(2*n,2*i+1,2*j+1)) += 0.25*v(IDX(n,i,j));
        }
    return res;
}

template class LinearInterpolation<1>;
template class LinearInterpolation<2>;