#pragma once

#include "Core/Function.h"
#include "Core/Norm.h"
#include "Core/matrix.h"
#include "Core/sparseMatrix.h"
#include "IntergridOp/IntergridOp.h"
#include "Boundary.h"
#include <iostream>
#include <fstream>
#include <iomanip>

template <int Dim>
class Solver{
protected:
    int n; // Grid Size
    ColVector cur, b;
    const IntergridOp<Dim> &restriction;
    const IntergridOp<Dim> &prolongation;
    std::string bondary;
    BType bonDetail;
    double eps;
    int maxiter;
    std::vector<SparseMatrix> Ah;
    int cycle; // cycle=1:V-Cycle  cycle=2:FMG-Cycle
    bool pure_Neumann;
    double neumannC;
    bool irregular;
    bool nineStencil;

    SparseMatrix getAh(const int &n) const;
    
    void initAh(){
        for(int i = n; i >= 2; i >>= 1)
            Ah.push_back(getAh(i));
    }

    ColVector VC(const int &d, const int &n, ColVector v, const ColVector &f);
    ColVector FMG(const int &d, const int &n, const ColVector &f);

public:
    Solver(const IntergridOp<Dim> &restriction, 
           const IntergridOp<Dim> &prolongation);
    
    void init(const int &_n, 
              const Function<Dim> &f, 
              const Function<Dim> &g, 
              const std::string &bon){
        init(_n, f, g, bon, BType());
    }

    void init(const int &_n, 
              const Function<Dim> &f, 
              const Function<Dim> &g, 
              const std::string &bon, 
              const BType &bonDetail);

    void setEps(const double &_eps) {eps = _eps;}
    void setMaxiter(const int &N) {maxiter = N;}
    void setCycle(const std::string &cy);
    void setIrregular(const bool &rhs) {irregular = rhs;}
    void useNineStencil() {nineStencil = true;}
    bool isPureNeumann() const {return pure_Neumann;}

    void solve();

    void output(std::ostream& fout) const;
    double calcNeumannC(const Function<Dim> &u);
    double checkError(const Function<Dim> &u, const Norm &norm) const;
};

double lowBon(const double &x); // Irregular down-bondary