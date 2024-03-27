#ifndef _V_CYCLE_H_
#define _V_CYCLE_H_

#include "matrix.h"
#include "sparseMatrix.h"
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

class Operator{
protected:
    int dim;
public:
    Operator(const int &_dim);
    virtual ColVector operator () (const ColVector &v) const = 0;
};

class injection : public Operator{
public:
    injection(const int &_dim);
    ColVector operator () (const ColVector &v) const;
};

class full_operator : public Operator{
public:
    full_operator(const int &_dim);
    ColVector operator () (const ColVector &v) const;
};

class linear_interpolation : public Operator{
public:
    linear_interpolation(const int &_dim);
    ColVector operator () (const ColVector &v) const;
};

class quadradic_interpolation : public Operator{
public:
    quadradic_interpolation(const int &_dim);
    ColVector operator () (const ColVector &v) const;
};


class Function{
public:
    virtual double operator () (const double &x) const = 0;
    virtual double operator () (const double &x, const double &y) const = 0;
    virtual double delta(const double &x, const double &y) const{return 0;}
};

class Function1D : public Function{
public:
    double operator () (const double &x, const double &y) const{
        return 0;
    }
};

class Function2D : public Function{
public:
    double operator () (const double &x) const{
        return 0;
    }
};

class BType{
private:
    std::map<std::string, std::string> types;
public:
    void setBondary(const std::string& bon, const std::string& typ);
    std::string operator () (const std::string& bon) const;
};

class Norm{
public:
    virtual double operator () (const std::vector<double> &vec) const = 0;
};

class rootSolver{
protected:
    int n; // Grid Size
    ColVector cur, b;
    const Operator & restriction;
    const Operator & prolongation;
    std::string bondary;
    BType bonDetail;
    double eps;
    int maxiter;
    int dim;
    std::vector<SparseMatrix> Ah;
    int cycle; // cycle=1:V-Cycle  cycle=2:FMG-Cycle
    bool pure_Neumann;
    double neumannC;
    bool irregular;
    bool nineStencil;

    virtual SparseMatrix getAh(const int &n) const = 0;
    void initAh();
    ColVector VC(const int &d, const int &n, ColVector v, const ColVector &f);
    ColVector FMG(const int &d, const int &n, const ColVector &f);

public:
    rootSolver(const Operator & restriction, const Operator & prolongation);
    void init(const int &_n, const Function &f, const Function &g, const std::string &bon);
    virtual void init(const int &_n, const Function &f, const Function &g, const std::string &bon, const BType &bonDetail) = 0;
    void setEps(const double &_eps);
    void setMaxiter(const int &N);
    void setCycle(const std::string &cy);
    void setIrregular(const bool &rhs);
    void useNineStencil();
    void solve();
    bool isPureNeumann() const;
    virtual void output(std::ostream& fout) const = 0;
    virtual double calcNeumannC(const Function &u) = 0;
    virtual double checkError(const Function &u, const Norm &norm) const = 0;
};

class Solver1D : public rootSolver{
private:
    SparseMatrix getAh(const int &n) const;
public:
    Solver1D(const Operator & restriction, const Operator & prolongation);
    void init(const int &_n, const Function &f, const Function &g, const std::string &bon, const BType &bonDetail);
    void output(std::ostream& fout) const;
    double calcNeumannC(const Function &u);
    double checkError(const Function &u, const Norm &norm) const;
};

class Solver2D : public rootSolver{
private:
    SparseMatrix getAh(const int &n) const;
public:
    Solver2D(const Operator & restriction, const Operator & prolongation);
    void init(const int &_n, const Function &f, const Function &g, const std::string &bon, const BType &bonDetail);
    void output(std::ostream& fout) const;
    double calcNeumannC(const Function &u);
    double checkError(const Function &u, const Norm &norm) const;
};

double lowBon(const double &x); // Irregular down-bondary

#endif