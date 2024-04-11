#include "QuadraticInterpolation.h"
#include "LinearInterpolation.h"

#include <cmath>

template<>
ColVector QuadradicInterpolation<1>::operator () (const ColVector &v) const{
    if(v.n <= 2){
        //点数少于2个无法二次插值
        LinearInterpolation<1> lin;
        return lin(v);
    }
    ColVector res(v.n*2-1);
    for(int i = 0; i < v.n; i++)
        res(2*i) = v(i);
    res(1) = 0.375*v(0) + 0.75*v(1) - 0.125*v(2);
    res(res.n-2) = 0.375*v(v.n-1) + 0.75*v(v.n-2) - 0.125*v(v.n-3);
    for(int i = 1; i < v.n-2; i++)
        res(2*i+1) = 9.0/16.0*(v(i)+v(i+1)) - 1.0/16.0*(v(i-1)+v(i+2));
    return res;
}

template<>
ColVector QuadradicInterpolation<2>::operator () (const ColVector &v) const{
    if(v.n <= 2){
        //点数少于2个无法二次插值
        LinearInterpolation<2> lin;
        return lin(v);
    }
    int n = (int)round(sqrt(v.n)) - 1;
    ColVector res( (2*n+1)*(2*n+1) );
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++)
            res(IDX(2*n,2*i,2*j)) = v(IDX(n,i,j));
    for(int i = 0; i <= n; i++){
        res(IDX(2*n,2*i,1)) = 0.375*v(IDX(n,i,0)) + 0.75*v(IDX(n,i,1)) - 0.125*v(IDX(n,i,2));
        res(IDX(2*n,2*i,2*n-1)) = 0.375*v(IDX(n,i,n)) + 0.75*v(IDX(n,i,n-1)) - 0.125*v(IDX(n,i,n-2));
        for(int j = 1; j < n-1; j++)
            res(IDX(2*n,2*i,2*j+1)) = 9.0/16.0*(v(IDX(n,i,j))+v(IDX(n,i,j+1))) - 1.0/16.0*(v(IDX(n,i,j-1))+v(IDX(n,i,j+2)));
    }
    for(int j = 0; j <= n; j++){
        res(IDX(2*n,1,2*j)) = 0.375*v(IDX(n,0,j)) + 0.75*v(IDX(n,1,j)) - 0.125*v(IDX(n,2,j));
        res(IDX(2*n,2*n-1,2*j)) = 0.375*v(IDX(n,n,j)) + 0.75*v(IDX(n,n-1,j)) - 0.125*v(IDX(n,n-2,j));
        for(int i = 1; i < n-1; i++)
            res(IDX(2*n,2*i+1,2*j)) = 9.0/16.0*(v(IDX(n,i,j))+v(IDX(n,i+1,j))) - 1.0/16.0*(v(IDX(n,i-1,j))+v(IDX(n,i+2,j)));
    }
    res(IDX(2*n,1,1)) = 0.5*(v(IDX(n,0,1))+v(IDX(n,1,0))) + 0.25*v(IDX(n,1,1)) - 0.125*(v(IDX(n,0,2))+v(IDX(n,2,0)));
    res(IDX(2*n,1,2*n-1)) = 0.5*(v(IDX(n,0,n-1))+v(IDX(n,1,n))) + 0.25*v(IDX(n,1,n-1)) - 0.125*(v(IDX(n,0,n-2))+v(IDX(n,2,n)));
    res(IDX(2*n,2*n-1,1)) = 0.5*(v(IDX(n,n-1,0))+v(IDX(n,n,1))) + 0.25*v(IDX(n,n-1,1)) - 0.125*(v(IDX(n,n-2,0))+v(IDX(n,n,2)));
    res(IDX(2*n,2*n-1,2*n-1)) = 0.5*(v(IDX(n,n,n-1))+v(IDX(n,n-1,n))) + 0.25*v(IDX(n,n-1,n-1)) - 0.125*(v(IDX(n,n,n-2))+v(IDX(n,n-2,n)));
    for(int i = 1; i < n-1; i++){
        res(IDX(2*n,1,2*i+1)) = 0.25*(v(IDX(n,0,i))+v(IDX(n,0,i+1))) + 0.375*(v(IDX(n,1,i))+v(IDX(n,1,i+1))) - 0.0625*(v(IDX(n,0,i-1))+v(IDX(n,0,i+2))+v(IDX(n,2,i))+v(IDX(n,2,i+1)));
        res(IDX(2*n,2*i+1,1)) = 0.25*(v(IDX(n,i,0))+v(IDX(n,i+1,0))) + 0.375*(v(IDX(n,i,1))+v(IDX(n,i+1,1))) - 0.0625*(v(IDX(n,i-1,0))+v(IDX(n,i+2,0))+v(IDX(n,i,2))+v(IDX(n,i+1,2)));
        res(IDX(2*n,2*n-1,2*i+1)) = 0.25*(v(IDX(n,n,i))+v(IDX(n,n,i+1))) + 0.375*(v(IDX(n,n-1,i))+v(IDX(n,n-1,i+1))) - 0.0625*(v(IDX(n,n,i-1))+v(IDX(n,n,i+2))+v(IDX(n,n-2,i))+v(IDX(n,n-2,i+1)));
        res(IDX(2*n,2*i+1,2*n-1)) = 0.25*(v(IDX(n,i,n))+v(IDX(n,i+1,n))) + 0.375*(v(IDX(n,i,n-1))+v(IDX(n,i+1,n-1))) - 0.0625*(v(IDX(n,i-1,n))+v(IDX(n,i+2,n))+v(IDX(n,i,n-2))+v(IDX(n,i+1,n-2)));
    }
    for(int i = 1; i < n-1; i++)
        for(int j = 1; j < n-1; j++){
            res(IDX(2*n,2*i+1,2*j+1)) = 5.0/16.0 * (v(IDX(n,i,j))+v(IDX(n,i+1,j))+v(IDX(n,i,j+1))+v(IDX(n,i+1,j+1)))
                                        -1.0/32.0 * (v(IDX(n,i,j-1))+v(IDX(n,i,j+2))+v(IDX(n,i-1,j))+v(IDX(n,i-1,j+1))+v(IDX(n,i+1,j-1))+v(IDX(n,i+1,j+2))+v(IDX(n,i+2,j))+v(IDX(n,i+2,j+1)));
        }
    return res;
}

template class QuadradicInterpolation<1>;
template class QuadradicInterpolation<2>;