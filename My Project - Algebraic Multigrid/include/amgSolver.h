#ifndef _AMG_SOLVER_H_
#define _AMG_SOLVER_H_

#include "sparseMatrix.h"
#include "matrix.h"
#include <vector>

class amgSolver{
private:
    std::vector<SparseMatrix> Ah, Ph, Rh;
    // Ah:算子  Ph:插值  Rh:限制

    ColVector VC(const int &d, ColVector x, const ColVector &rhs);
    void generateGrid(const int &d);
    SparseMatrix getInterpolator(const SparseMatrix &A);
    std::vector<int> getCorsetPoints(const SparseMatrix &A);

public:
    void generateGrid(const SparseMatrix &A);
    ColVector solve(const ColVector &rhs, const int &maxIter, const double &eps);
};

#endif