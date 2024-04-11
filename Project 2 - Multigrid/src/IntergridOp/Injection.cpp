#include "Injection.h"

template<>
ColVector Injection<1>::operator () (const ColVector &v) const{
    ColVector res(v.n/2+1);
    for(int i = 0; i < res.n; i++){
        res(i) = v(2*i);
    }
    return res;
}

template<>
ColVector Injection<2>::operator () (const ColVector &v) const{
    int n = (int)round(sqrt(v.n)) - 1;
    ColVector res( (n/2+1)*(n/2+1) );
    for(int i = 0; i <= n/2; i++)
        for(int j = 0; j <= n/2; j++)
            res(IDX(n/2,i,j)) = v(IDX(n,2*i,2*j));
    return res;
}

template class Injection<1>;
template class Injection<2>;