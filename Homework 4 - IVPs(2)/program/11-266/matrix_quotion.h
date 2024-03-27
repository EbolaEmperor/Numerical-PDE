/***************************************************************
 *
 * 这是一个利用高精度分数运算实现的精确矩阵运算库
 * 运算效率较低，但能保证绝对精确的运算，可用于理论研究
 * 依赖于：fraction.h, bigint.h
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#ifndef _MATRIX_QUOTION_H_
#define _MATRIX_QUOTION_H_

#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "fraction.h"
#include "extendQuotion.h"

template<class Type> class Matrix{
private:
    std::vector<Type> a;
public:
    int n, m;
    Matrix(){
        n = m = 0;
        a.clear();
    }
    Matrix(const int &_n, const int &_m){
        n = _n;
        m = _m;
        a.resize(n*m);
        for(int i = 0; i < n*m; i++)
            a[i] = Type(0);
    }
    Matrix(const Matrix &A){
        n = A.n;
        m = A.m;
        a.resize(n*m);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                a[i*m+j] = A.element(i,j);
    }
    ~Matrix(){
        a.clear();
    }

    bool iszero() const{
        for(int i = 0; i < n*m; i++)
            if(a[i] != 0) return false;
        return true; 
    }
    Matrix & operator = (const Matrix & rhs){
        Matrix copy(rhs);
        std::swap(*this, copy);
        return *this;
    }
    Matrix & operator = (Matrix && rhs){
        std::swap(a, rhs.a);
        n = rhs.n;
        m = rhs.m;
        return *this;
    }

    // 矩阵基本运算
    Type & element(const int &i, const int &j){
        return a[i*m+j];
    }
    const Type & element(const int &i, const int &j) const{
        return a[i*m+j];
    }
    Type & operator () (const int &i, const int &j){
        return element(i, j);
    }
    const Type & operator () (const int &i, const int &j) const{
        return element(i, j);
    }

    friend std::ostream& operator << (std::ostream& out, const Matrix &A){
        for(int i = 0; i < A.n; i++)
        {
            out << "[ " << A.element(i,0);
            for(int j = 1; j < A.m; j++)
                out << ", " << A.element(i,j);
            out << " ]" << std::endl;
        }
        return out;
    }

    Matrix operator * (const Matrix &rhs) const{
        Matrix C(n, rhs.m);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < rhs.m; j++)
                for(int k = 0; k < m; k++)
                    C(i,j) += element(i,k) * rhs(k,j);
        return C;
    }

    Matrix operator - (const Matrix &rhs) const{
        Matrix C(n, m);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                    C(i,j) = element(i,j) - rhs(i,j);
        return C;
    }

    void swaprow(const int &r1, const int &r2){
        if(r1<0 || r1>=n || r2<0 || r2>=n){
            std::cerr << "Matrix Error! Swaprow ouof range!" << std::endl;
            exit(-1);
        }
        for(int j = 0; j < m; j++)
            std::swap(a[r1*m+j], a[r2*m+j]);
    }

    void swapcol(const int &r1, const int &r2){
        if(r1<0 || r1>=m || r2<0 || r2>=m){
            std::cerr << "Matrix Error! Swapcol ouof range!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < n; i++)
            std::swap(a[i*m+r1], a[i*m+r2]);
    }

    // 矩阵常用算法
    // friend Matrix solve(Matrix A, Matrix b){
    //     if(A.m!=A.n || A.n!=b.n || A.m==0){
    //         std::cerr << "Matrix Error! The method solve() cannot solve an ill-posed equation!" << std::endl;
    //         return Matrix();
    //     }
    //     int n = A.n;
    //     Matrix x(n,1);
    //     for(int i = 0; i < n; i++){
    //         int p = i;
    //         while(p<n && A.element(p,i)==0) p++;
    //         if(p==n){
    //             std::cerr << "Matrix Error! The method solve() cannot solve an singular equation!" << std::endl;
    //             return Matrix();
    //         }
    //         if(p!=i) A.swaprow(i,p), b.swaprow(i,p);
    //         for(int j = 0; j < n; j++){
    //             if(i==j) continue;
    //             Type coef = A.element(j,i)/A.element(i,i);
    //             for(int k = i; k < n; k++)
    //                 A.element(j,k) -= A.element(i,k)*coef;
    //             b.element(j,0) -= b.element(i,0)*coef;
    //         }
    //     }
    //     for(int i = 0; i < n; i++)
    //         x.element(i,0) = b.element(i,0)/A.element(i,i);
    //     return x;
    // }

    friend Type det(const Matrix &A){
        const int n = A.n;
        int *perm = new int[n];
        for(int i = 0; i < n; i++)
            perm[i] = i;
        Type res(0);
        do{
            Type mul = Type::one;
            for(int i = 0; i < n; i++)
                for(int j = i+1; j < n; j++)
                    if(perm[i] > perm[j]) mul = -mul;
            for(int i = 0; i < n; i++)
                mul = mul * A(i,perm[i]);
            res = res + mul;
        } while(std::next_permutation(perm, perm+n));
        return res;
    }
};

#endif
