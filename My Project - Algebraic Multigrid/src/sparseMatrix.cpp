#include "sparseMatrix.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <unordered_map>
using namespace std;

SparseMatrix::SparseMatrix(){
    n = m = size = 0;
    row_index = nullptr;
    elements = nullptr;
}

SparseMatrix::SparseMatrix(const int &_n, const int &_m){
    n = _n;
    m = _m;
    size = 0;
    row_index = nullptr;
    elements = nullptr;
}

int SparseMatrix::rowSize() const{
    return n;
}

int SparseMatrix::colSize() const{
    return m;
}

double SparseMatrix::density() const{
    return (double)size/((double)n*m);
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

SparseMatrix::SparseMatrix(const Matrix & rhs){
    n = rhs.n;
    m = rhs.m;
    size = 0;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            if(rhs(i,j)) size++;
    elements = new SparseElement[size];
    row_index = new int[n+1];
    int idx = 0;
    for(int i = 0; i < n; i++){
        row_index[i] = idx;
        for(int j = 0; j < m; j++)
            if(rhs(i,j)) elements[idx++] = SparseElement(j, rhs(i,j));
    }
    row_index[n] = idx;
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
    ColVector res(n);
    for(int i = 0; i < n; i++){
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
    RowVector res(A.m);
    for(int i = 0; i < A.m; i++){
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

ColVector SparseMatrix::GaussSeidel(ColVector x, const ColVector &b) const{
    for(int i = 0; i < n; i++){
        double coef = 0, sum = 0;
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            if(elements[c].j!=i) sum += elements[c].value * x(elements[c].j);
            else coef = elements[c].value;
        x(i) = (b(i) - sum) / coef;
    }
    return x;
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

SparseMatrix SparseMatrix::T() const{
    SparseMatrix res(m, n);
    vector<SparseElement> spe[m];
    for(int i = 0; i < n; i++)
        for(int c = row_index[i]; c < row_index[i+1]; c++)
            spe[elements[c].j].push_back(SparseElement(i, elements[c].value));
    res.size = size;
    res.elements = new SparseElement[size];
    res.row_index = new int[m+1];
    int idx = 0;
    for(int i = 0; i < m; i++){
        res.row_index[i] = idx;
        for(auto& el : spe[i])
            res.elements[idx++] = el;
    }
    res.row_index[m] = idx;
    return res;
}

vector<int> SparseMatrix::nonzeroIndexInRow(const int &r) const{
    vector<int> p;
    for(int c = row_index[r]; c < row_index[r+1]; c++)
        p.push_back(elements[c].j);
    return p;
}

vector<double> SparseMatrix::nonzeroValueInRow(const int &r) const{
    vector<double> p;
    for(int c = row_index[r]; c < row_index[r+1]; c++)
        p.push_back(elements[c].value);
    return p;
}

double SparseMatrix::operator () (const int &i, const int &j) const{
    for(int c = row_index[i]; c < row_index[i+1]; c++){
        if(elements[c].j == j){
            return elements[c].value;
        } else if(elements[c].j > j){
            break;
        }
    }
    return 0;
}

SparseMatrix SparseMatrix::operator * (const SparseMatrix &rhs) const{
    unordered_map<long long, double> f;
    for(int i = 0; i < n; i++)
        for(int c = row_index[i]; c < row_index[i+1]; c++){
            int k = elements[c].j;
            for(int s = rhs.row_index[k]; s < rhs.row_index[k+1]; s++){
                f[ 1ll*i*rhs.m + rhs.elements[s].j ] += elements[c].value * rhs.elements[s].value;
            }
        }
    vector<Triple> elem;
    for(auto p : f){
    	if(fabs(p.second)<1e-16) continue;
        int r = p.first / rhs.m;
        int c = p.first - 1ll*rhs.m*r;
        elem.push_back(Triple(r,c,p.second));
    }
    return SparseMatrix(n, rhs.m, elem);
}
