/***************************************************************
 *
 * 这是一个矩阵运算库，为了方便以后设计算法更加简洁，特编写以用
 * 版本号：v1.0.2
 * 
 * copyright © 2022 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cstring>
#include <cmath>

const double mat_eps = 1e-12;

class Matrix;
class ColVector;
class RowVector;

class Matrix{
private:
    double *a;
public:
    int n, m;
    Matrix();
    Matrix(const int &_n);
    Matrix(const int &_n, const int &_m);
    Matrix(const Matrix &A);
    Matrix(const double *p, const int &_n);
    Matrix(const int &_n, const int &_m, const double *p);
    ~Matrix();
    bool empty() const;

    Matrix & operator = (const Matrix & rhs);
    Matrix & operator = (Matrix && rhs);

    const double operator () (const int &r, const int &c) const;
    double & operator () (const int &r, const int &c);
    const double element(const int &r, const int &c) const;
    double & element(const int &r, const int &c);

    friend Matrix diag(const Matrix &A);

    void setSubmatrix(const int &r, const int &c, const Matrix &rhs);
    Matrix getSubmatrix(const int &r1, const int &r2, const int &c1, const int &c2) const;
    Matrix reshape(const int &_n, const int &_m) const;

    RowVector getRow(const int &r) const;
    ColVector getCol(const int &c) const;

    Matrix operator + (const Matrix &B) const;
    Matrix operator - () const;
    Matrix operator - (const Matrix &B) const;
    Matrix operator * (const Matrix &B) const;
    Matrix operator / (const double &p) const;
    Matrix T() const;

    double vecnorm(const double &p) const;
    void swaprow(const int &r1, const int &r2);
    void swapcol(const int &r1, const int &r2);

    friend Matrix solve(Matrix A, Matrix b);
    ColVector solve(const ColVector &b) const;
    double det() const;
    Matrix inv() const;
    Matrix rref() const;
    void FGdecompose(Matrix &F, Matrix &G) const;
    Matrix pinv() const;
    double sqrsum() const;
};

class RowVector: public Matrix{
public:
    RowVector(): Matrix() {};
    RowVector(const int &n): Matrix(1,n) {};
    RowVector(const int &n, const double *p): Matrix(1,n,p) {};
    RowVector(const Matrix &rhs);
    int size() const;
    const double operator ()(const int &x) const;
    double & operator () (const int &x);
    RowVector operator + (const RowVector &rhs) const;
    RowVector operator - (const RowVector &rhs) const;
    RowVector operator - () const;
    ColVector T() const;
};

class ColVector: public Matrix{
public:
    ColVector(): Matrix() {};
    ColVector(const int &n): Matrix(n,1) {};
    ColVector(const int &n, const double *p): Matrix(n,1,p) {};
    ColVector(const Matrix &rhs);
    int size() const;
    const double operator ()(const int &x) const;
    double & operator () (const int &x);
    ColVector operator + (const ColVector &rhs) const;
    ColVector operator - (const ColVector &rhs) const;
    ColVector operator - () const;
    RowVector T() const;
};

Matrix hilbert(const int &n);
Matrix zeros(const int &n, const int &m);
Matrix ones(const int &n, const int &m);
Matrix eye(const int &n);
double value(const Matrix &A);
ColVector zeroCol(const int &n);
RowVector zeroRow(const int &n);
int sgn(const double &x);

//----------------------Matrix相关函数---------------------------
Matrix operator * (const double &k, const Matrix &x);
Matrix abs(const Matrix &A);
double max(const Matrix &A);
double sum(const Matrix &A);
std::istream& operator >> (std::istream& in, Matrix &A);
std::ostream& operator << (std::ostream& out, const Matrix &A);
Matrix dotdiv(const Matrix &a, const Matrix &b);
Matrix solveLowerTriangular(const Matrix &A, const Matrix &b);
Matrix solveUpperTriangular(const Matrix &A, const Matrix &b);
ColVector CG_solve(const Matrix &A, const ColVector &b);
ColVector CG_solve(const Matrix &A, const ColVector &b, const double err);
ColVector CG_solve(const Matrix &A, const ColVector &b, const double err, ColVector x);
double det(const Matrix &A);
Matrix inv(const Matrix &A);
Matrix choleskyImproved(const Matrix &A);
Matrix solveByLDL(const Matrix &A, const Matrix &b);
Matrix gillMurray(Matrix A);
Matrix solveByLDL_GM(const Matrix &A, const Matrix &b);
Matrix pinv(const Matrix &A);
Matrix mergeCol(const Matrix &A, const Matrix &B);
Matrix mergeRow(const Matrix &A, const Matrix &B);
Matrix min(const Matrix &A, const Matrix &B);
Matrix max(const Matrix &A, const Matrix &B);
double vecnorm(const Matrix &A, const double &p);
double vecnorm(const Matrix &A);

//----------------------Row/ColVector相关函数----------------------
RowVector operator * (const double &k, const RowVector &x);
ColVector operator * (const double &k, const ColVector &x);
ColVector operator * (const Matrix &A, const ColVector &x);
RowVector operator * (const RowVector &x, const Matrix &A);
double operator * (const RowVector &r, const ColVector &c);

inline int IDX(const int &n, const int &i, const int &j){
    return (n + 1) * i + j;
}

inline ColVector zeroGridCol(const int &n, const int &dim){
    return zeros( (int)round(pow(n+1,dim)), 1 );
}

#endif