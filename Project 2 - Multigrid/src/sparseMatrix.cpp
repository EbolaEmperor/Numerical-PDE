#include "sparseMatrix.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>
using namespace std;

SparseMatrix::SparseMatrix(){
    n = m = size = 0;
    row_index = nullptr;
    elements = nullptr;
}

SparseMatrix::SparseMatrix(const int &_n, const int &_m, vector<Triple> &ele){
    n = _n;
    m = _m;
    row_index = new int[n+1];
    sort(ele.begin(), ele.end());
    int num = 0;
    for(int i = 0; i < ele.size(); i++){
        num++;
        if(i && ele[i].i==ele[i-1].i && ele[i].j==ele[i-1].j)
            num--;
    }
    elements = new SparseElement[num];
    size = num;
    num = 0;
    int row = 0;
    row_index[0] = 0;
    for(int i = 0; i < ele.size(); i++){
        if(i && ele[i].i==ele[i-1].i && ele[i].j==ele[i-1].j){
            elements[num].value += ele[i].value;
        } else {
            while(ele[i].i > row)
                row_index[++row] = num;
            elements[num] = SparseElement(ele[i].j, ele[i].value);
            num++;
        }
    }
    while(row < n)
        row_index[++row] = size;
}

SparseMatrix::SparseMatrix(const SparseMatrix & rhs){
    n = rhs.n;
    m = rhs.m;
    size = rhs.size;
    row_index = new int[n+1];
    memcpy(row_index, rhs.row_index, sizeof(int)*(n+1));
    elements = new SparseElement[size];
    memcpy(elements, rhs.elements, sizeof(SparseElement)*size);
}

SparseMatrix::~SparseMatrix(){
    clear();
}

void SparseMatrix::clear(){
    n = m = size = 0;
    delete [] row_index;
    delete [] elements;
    row_index = nullptr;
    elements = nullptr;
}

SparseMatrix SparseMatrix::operator + (const SparseMatrix &rhs) const{
    if(n!=rhs.n || m!=rhs.m){
        cerr << "[Error] Cannot use operator + at matrixs of distinct size!" << endl;
        exit(-1);
    }
    vector<Triple> vec;
    for(int i = 0; i < n; i++)
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            vec.push_back(Triple(i, elements[j].j, elements[j].value));
    for(int i = 0; i < n; i++)
        for(int j = rhs.row_index[i]; j < rhs.row_index[i+1]; j++)
            vec.push_back(Triple(i, rhs.elements[j].j, rhs.elements[j].value));
    return SparseMatrix(n, m, vec);
}

SparseMatrix SparseMatrix::operator - (const SparseMatrix &rhs) const{
    if(n!=rhs.n || m!=rhs.m){
        cerr << "[Error] Cannot use operator + at matrixs of distinct size!" << endl;
        exit(-1);
    }
    vector<Triple> vec;
    for(int i = 0; i < n; i++)
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            vec.push_back(Triple(i, elements[j].j, elements[j].value));
    for(int i = 0; i < n; i++)
        for(int j = rhs.row_index[i]; j < rhs.row_index[i+1]; j++)
            vec.push_back(Triple(i, rhs.elements[j].j, -rhs.elements[j].value));
    return SparseMatrix(n, m, vec);
}

ColVector SparseMatrix::operator * (const ColVector & rhs) const{
    if(m!=rhs.n){
        cerr << "[Error] The columns of SparseMatrix does not coincide the rows of ColVector!" << endl;
        exit(-1);
    }
    ColVector res(rhs.n);
    for(int i = 0; i < rhs.n; i++){
        for(int j = row_index[i]; j < row_index[i+1]; j++)
            res(i) += elements[j].value * rhs(elements[j].j);
    }
    return res;
}

RowVector operator * (const RowVector & lhs, const SparseMatrix &A){
    if(A.n!=lhs.m){
        cerr << "[Error] The rows of SparseMatrix does not coincide the columns of RowVector!" << endl;
        exit(-1);
    }
    RowVector res(lhs.n);
    for(int i = 0; i < lhs.n; i++){
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            res(A.elements[j].j) += A.elements[j].value * lhs(i);
    }
    return res;
}

ColVector SparseMatrix::wJacobi(const ColVector &x, const ColVector &b, const double &w) const{
    ColVector x1 = b;
    for(int i = 0; i < n; i++){
        double coef = 0;
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            if(elements[c].j!=i) x1(i) -= elements[c].value * x(elements[c].j);
            else coef = elements[c].value;
        x1(i) /= coef;
    }
    return (1-w)*x + w*x1;
}

ostream & operator << (std::ostream & out, const SparseMatrix &A){
    out << "shape: " << A.n << " * " << A.m << endl;
    out << "non-zero elements:" << endl;
    for(int i = 0; i < A.n; i++)
        for(int j = A.row_index[i]; j < A.row_index[i+1]; j++)
            out << "(" << i << ", " << A.elements[j].j << ", " << A.elements[j].value << ")"<< std::endl;
    out << "row_index:" << endl;
    for(int i = 0; i <= A.n; i++)
        out << A.row_index[i] << ", ";
    out << endl;
    return out;
}

Matrix SparseMatrix::toDense() const{
    Matrix A(n,m);
    for(int i = 0; i < n; i++)
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            A(i, elements[c].j) = elements[c].value;
    return A;
}

ColVector SparseMatrix::LUsolve(const ColVector &b) const{
    Matrix A = toDense();
    return A.solve(b);
}