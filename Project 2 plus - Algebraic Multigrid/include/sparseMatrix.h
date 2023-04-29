#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include "matrix.h"
#include <vector>

struct SparseElement{
    int j;
    double value;
    SparseElement(const int &j=0, const double &value=0):
        j(j), value(value){}
};

struct Triple{
    int i, j;
    double value;
    Triple(): i(0), j(0), value(0.0) {}
    Triple(const int &i, const int &j, const double &value):
        i(i), j(j), value(value){}
    ~Triple(){}
    bool operator < (const Triple &b) const{
        return i==b.i && j<b.j || i<b.i;
    }
};

class SparseMatrix{
private:
    int n, m, size;
    int * row_index;
    SparseElement *elements;

public:
    SparseMatrix();
    SparseMatrix(const int &_n, const int &_m);
    SparseMatrix(const int &_n, const int &_m, std::vector<Triple> &ele);
    SparseMatrix(const SparseMatrix & rhs);
    SparseMatrix(const Matrix & rhs);
    ~SparseMatrix();
    void clear();
    int rowSize() const;
    int colSize() const;
    double density() const;
    SparseMatrix T() const;
    SparseMatrix operator * (const SparseMatrix &rhs) const;
    std::vector<int> nonzeroIndexInRow(const int &r) const;
    std::vector<double> nonzeroValueInRow(const int &r) const;
    double operator () (const int &i, const int &j) const;
    ColVector operator * (const ColVector &rhs) const;
    friend RowVector operator * (const RowVector &lhs, const SparseMatrix &A);
    SparseMatrix operator + (const SparseMatrix &rhs) const;
    SparseMatrix operator - (const SparseMatrix &rhs) const;
    ColVector wJacobi(const ColVector & x, const ColVector & b, const double &w) const;
    ColVector GaussSeidel(ColVector x, const ColVector & b) const;
    friend std::ostream & operator << (std::ostream & out, const SparseMatrix &A);
    Matrix toDense() const;
    ColVector LUsolve(const ColVector &b) const;
};

#endif