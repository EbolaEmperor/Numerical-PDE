#include "solver.h"
#include "matrix.h"
#include "sparseMatrix.h"
#include <iostream>

//-----------------------------BType------------------------------------

void BType::setBondary(const std::string& bon, const std::string& typ){
    if(typ != "Dirichlet" && typ != "Neumann"){
        std::cerr << "Unrecignized bondary '" << typ << "'" << std::endl;
        exit(-1);
    }
    types[bon] = typ;
}

std::string BType:: operator () (const std::string& bon) const{
    if(!types.count(bon)){
        std::cerr << "Undefined bondary '" << bon << "'" << std::endl;
        exit(-1);
    }
    return types.find(bon)->second;
}

//--------------------------------Basic Tools---------------------------

const int VC_V1 = 3;
const int VC_V2 = 3;
const int VC_MAX_ITER = 100;
const double VC_W = 2.0/3.0;
const double VC_PI = acos(-1);
const double VC_EPS = 1e-15;

double lowBon(const double &x){
    return sin(VC_PI*x)/16;
}

ColVector zeroGridCol(const int &n, const int &dim){
    return zeros( (int)round(pow(n+1,dim)), 1 );
}

void checkDim(const int &dim){
    if(dim!=1 && dim!=2){
        std::cerr << "[Error] dimension parameter can only be 1 or 2." << std::endl;
        exit(-1);
    }
}

int IDX(const int &n, const int &i, const int &j){
    return (n + 1) * i + j;
}

struct specialPoint{
    double x, y;
    specialPoint(): x(0.0), y(0.0) {}
    specialPoint(const int &n, const int &i, const int &j){
        double h = 1.0/n;
        double _x = i*h, _y = j*h;
        x = _x;
        y = lowBon(x) + _y * (1.0-lowBon(x));
    }
};

//--------------------------------Operators-----------------------------

Operator::Operator(const int &_dim){
    dim = _dim;
}

//--------------------------------Injection-----------------------------

injection::injection(const int &_dim):Operator(_dim){}

ColVector injection::operator () (const ColVector &v) const{
    if(dim==1){
        ColVector res(v.n/2+1);
        for(int i = 0; i < res.n; i++){
            res(i) = v(2*i);
        }
        return res;
    } else {
        int n = (int)round(sqrt(v.n)) - 1;
        ColVector res( (n/2+1)*(n/2+1) );
        for(int i = 0; i <= n/2; i++)
            for(int j = 0; j <= n/2; j++)
                res(IDX(n/2,i,j)) = v(IDX(n,2*i,2*j));
        return res;
    }
}

//------------------------------Full Operator--------------------------------

full_operator::full_operator(const int &_dim):Operator(_dim){}

ColVector full_operator::operator () (const ColVector &v) const{
    if(dim==1){
        ColVector res(v.n/2+1);
        for(int i = 1; i < res.n-1; i++){
            res(i) = 0.25*v(2*i-1) + 0.5*v(2*i) + 0.25*v(2*i+1);
        }
        res(0) = v(0);
        res(v.n/2) = v(v.n-1);
        return res;
    } else {
        int n = (int)round(sqrt(v.n)) - 1;
        ColVector res( (n/2+1)*(n/2+1) );
        res(IDX(n/2,0,0)) = v(IDX(n,0,0));
        res(IDX(n/2,0,n/2)) = v(IDX(n,0,n));
        res(IDX(n/2,n/2,0)) = v(IDX(n,n,0));
        res(IDX(n/2,n/2,n/2)) = v(IDX(n,n,n));
        for(int i = 1; i < n/2; i++){
            res(IDX(n/2,i,0)) = 0.25*v(IDX(n,2*i-1,0)) + 0.5*v(IDX(n,2*i,0)) + 0.25*v(IDX(n,2*i+1,0));
            res(IDX(n/2,i,n/2)) = 0.25*v(IDX(n,2*i-1,n)) + 0.5*v(IDX(n,2*i,n)) + 0.25*v(IDX(n,2*i+1,n));
            res(IDX(n/2,0,i)) = 0.25*v(IDX(n,0,2*i-1)) + 0.5*v(IDX(n,0,2*i)) + 0.25*v(IDX(n,0,2*i+1));
            res(IDX(n/2,n/2,i)) = 0.25*v(IDX(n,n,2*i-1)) + 0.5*v(IDX(n,n,2*i)) + 0.25*v(IDX(n,n,2*i+1));
        }
        for(int i = 1; i < n/2; i++)
            for(int j = 1; j < n/2; j++){
                res(IDX(n/2,i,j)) = 0.0625 * (v(IDX(n,2*i-1,2*j-1)) + v(IDX(n,2*i-1,2*j+1)) + v(IDX(n,2*i+1,2*j-1)) + v(IDX(n,2*i+1,2*j+1)))
                                   + 0.125 * (v(IDX(n,2*i-1,2*j)) + v(IDX(n,2*i+1,2*j)) + v(IDX(n,2*i,2*j-1)) + v(IDX(n,2*i,2*j+1)))
                                   +  0.25 * v(IDX(n,2*i,2*j));
            }
        return res;
    }
}

//---------------------------Linear Interpolation-----------------------------

linear_interpolation::linear_interpolation(const int &_dim):Operator(_dim){}

ColVector linear_interpolation::operator () (const ColVector &v) const{
    if(dim==1){
        ColVector res(v.n*2-1);
        for(int i = 0; i < v.n; i++){
            res(2*i) = v(i);
            if(i) res(2*i-1) += 0.5*v(i);
            if(i<v.n-1) res(2*i+1) += 0.5*v(i);
        }
        return res;
    } else {
        int n = (int)round(sqrt(v.n)) - 1;
        ColVector res( (2*n+1)*(2*n+1) );
        for(int i = 0; i <= n; i++)
            for(int j = 0; j <= n; j++){
                res(IDX(2*n,2*i,2*j)) = v(IDX(n,i,j));
                if(i) res(IDX(2*n,2*i-1,2*j)) += 0.5*v(IDX(n,i,j));
                if(i<n) res(IDX(2*n,2*i+1,2*j)) += 0.5*v(IDX(n,i,j));
                if(j) res(IDX(2*n,2*i,2*j-1)) += 0.5*v(IDX(n,i,j));
                if(j<n) res(IDX(2*n,2*i,2*j+1)) += 0.5*v(IDX(n,i,j));
                if(i && j) res(IDX(2*n,2*i-1,2*j-1)) += 0.25*v(IDX(n,i,j));
                if(i && j<n) res(IDX(2*n,2*i-1,2*j+1)) += 0.25*v(IDX(n,i,j));
                if(i<n && j) res(IDX(2*n,2*i+1,2*j-1)) += 0.25*v(IDX(n,i,j));
                if(i<n && j<n) res(IDX(2*n,2*i+1,2*j+1)) += 0.25*v(IDX(n,i,j));
            }
        return res;
    }
}

//---------------------------Quadradic Interpolation-----------------------------

quadradic_interpolation::quadradic_interpolation(const int &_dim):Operator(_dim){}

ColVector quadradic_interpolation::operator () (const ColVector &v) const{
    if(v.n <= 2){
        //点数少于2个无法二次插值
        linear_interpolation lin(dim);
        return lin(v);
    }
    if(dim==1){
        ColVector res(v.n*2-1);
        for(int i = 0; i < v.n; i++)
            res(2*i) = v(i);
        res(1) = 0.375*v(0) + 0.75*v(1) - 0.125*v(2);
        res(res.n-2) = 0.375*v(v.n-1) + 0.75*v(v.n-2) - 0.125*v(v.n-3);
        for(int i = 1; i < v.n-2; i++)
            res(2*i+1) = 9.0/16.0*(v(i)+v(i+1)) - 1.0/16.0*(v(i-1)+v(i+2));
        return res;
    } else {
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
}

//--------------------------------rootSolver--------------------------------

rootSolver::rootSolver(const Operator & restriction, const Operator & prolongation)
    : restriction(restriction), prolongation(prolongation){
    eps = 1e-16;
    maxiter = VC_MAX_ITER;
    cycle = 0;
    pure_Neumann = false;
    irregular = false;
    nineStencil = false;
}

void rootSolver::setCycle(const std::string &cy){
    if(cy=="V") cycle = 1;
    else if(cy=="FMG") cycle = 2;
    else{
        std::cerr << "[Error] Unrecignized cycle type '" << cy << "'" << std::endl;
        exit(-1);
    }
}

void rootSolver::setEps(const double &_eps){
    eps = _eps;
}

void rootSolver::setMaxiter(const int &N){
    maxiter = N;
}

bool rootSolver::isPureNeumann() const{
    return pure_Neumann;
}

void rootSolver::init(const int &_n, const Function &f, const Function &g, const std::string &bon){
    init(_n, f, g, bon, BType());
}

ColVector rootSolver::VC(const int &d, const int &n, ColVector v, const ColVector &f){
    const SparseMatrix & A = Ah[d];
    for(int i = 0; i < VC_V1; i++)
        v = A.wJacobi(v, f, VC_W);
    if(n > 2){
        v = v + prolongation( VC(d+1, n/2, zeroGridCol(n/2,dim), restriction(f-A*v)) );
    }
    for(int i = 0; i < VC_V2; i++)
        v = A.wJacobi(v, f, VC_W);
    return v;
}

ColVector rootSolver::FMG(const int &d, const int &n, const ColVector &f){
    ColVector v;
    if(n > 2){
        v = prolongation( FMG(d+1, n/2, restriction(f)) );
    } else {
        v = zeroGridCol(n,dim);
    }
    return VC(d, n, v, f);
}

void rootSolver::initAh(){
    // 初始化每一层的Ah，以减少VC过程中重复申请大量内存
    for(int i = n; i >= 2; i >>= 1){
        Ah.push_back(getAh(i));
    }
}

void rootSolver::solve(){
    std::cout << "Preparing coefficient matrixs..." << std::endl;
    initAh();
    const SparseMatrix &A = Ah[0];
    ColVector upd;
    int T = 0;
    double lastres = 1e100, curres;
    do{
        T++;
        std::cout << "Iteration " << T << "...";
        if(cycle==1)
            upd = VC(0, n, zeroGridCol(n,dim), b-A*cur);
        else if(cycle==2)
            upd = FMG(0, n, b-A*cur);
        else{
            std::cerr << "[Error] Missing definition of cycle-type before solve!" << std::endl;
            exit(-1);
        }
        curres = max(abs(upd));
        if(curres > lastres){
            std::cout << "   Terminated." << std::endl;
            break;
        }
        std::cout << "   Updation magnitude: " << curres << std::endl;
        cur = cur + upd;
        lastres = curres;
    }while(T < maxiter && curres > eps);
    Ah.clear();
}

void rootSolver::setIrregular(const bool &rhs){
    irregular = rhs;
}

void rootSolver::useNineStencil(){
    nineStencil = true;
}

//--------------------------------Solver 1D--------------------------------

Solver1D::Solver1D(const Operator & restriction, const Operator & prolongation)
    : rootSolver(restriction, prolongation){
    dim = 1;
}

void Solver1D::init(const int &_n, const Function &f, const Function &g, const std::string &bon, const BType &_bonDetail){
    n = _n;
    bondary = bon;
    bonDetail = _bonDetail;
    cur = zeroGridCol(n,1);
    b = zeroGridCol(n,1);
    double h = 1.0/n;
    for(int i = 1; i < n; i++){
        b(i) = f(i*h);
    }
    b(0) = g(0);
    b(n) = g(1);
    if(bon=="Neumann" || bon=="mixed" && bonDetail("x=0")=="Neumann" && bonDetail("x=1")=="Neumann")
        pure_Neumann = true;
}

SparseMatrix Solver1D::getAh(const int &n) const{
    double h = 1.0/n;
    std::vector<Triple> ele;
    for(int i = 1; i < n; i++){
        ele.emplace_back(Triple(i, i, 2.0/(h*h)));
        ele.emplace_back(Triple(i, i-1, -1.0/(h*h)));
        ele.emplace_back(Triple(i, i+1, -1.0/(h*h)));
    }
    if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("x=0") == "Dirichlet"){
        ele.emplace_back(Triple(0, 0, 1));
    } else {
        ele.emplace_back(Triple(0, 0, -1.5/h));
        ele.emplace_back(Triple(0, 1, 2.0/h));
        ele.emplace_back(Triple(0, 2, -0.5/h));
    }
    if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("x=1") == "Dirichlet"){
        ele.emplace_back(Triple(n, n, 1));
    } else {
        ele.emplace_back(Triple(n, n-2, 0.5/h));
        ele.emplace_back(Triple(n, n-1, -2.0/h));
        ele.emplace_back(Triple(n, n, 1.5/h));
    }
    return SparseMatrix(n+1, n+1, ele);
}

void Solver1D::output(std::ostream &out) const{
    out << std::fixed << std::setprecision(16);
    for(int i = 0; i <= n; i++)
        out << 1.0*i/n << " " << cur(i) << "\n";
}

double Solver1D::calcNeumannC(const Function &u){
    double l = 1e100, r = -1e100;
    double h = 1.0/n;
    for(int i = 0; i <= n; i++){
        l = std::min(l, u(i*h) - cur(i));
        r = std::max(r, u(i*h) - cur(i));
    }
    return neumannC = (l+r)/2;
}

double Solver1D::checkError(const Function &u, const Norm &norm) const{
    std::vector<double> vals;
    double h = 1.0/n;
    for(int i = 0; i <= n; i++)
        vals.emplace_back( cur(i) + neumannC - u(i*h) );
    return norm(vals);
}

//-------------------Special Tools for Irregular Grid----------------------

ColVector normalVector(const double &x){
    ColVector n(2);
    n(0) = VC_PI/16*cos(VC_PI*x);
    n(1) = -1;
    return n / vecnorm(n);
}

Matrix specialGenA(const std::vector<specialPoint> &sp){
    Matrix A(6,6);
    for(int i = 0; i < 6; i++){
        A(0,i) = 1;
        A(1,i) = sp[i].x - sp[0].x;
        A(2,i) = sp[i].y - sp[0].y;
        A(3,i) = 0.5 * (sp[i].x - sp[0].x) * (sp[i].x - sp[0].x);
        A(4,i) = 0.5 * (sp[i].y - sp[0].y) * (sp[i].y - sp[0].y);
        A(5,i) = (sp[i].x - sp[0].x) * (sp[i].y - sp[0].y);
    }
    return A;
}

ColVector specialSolveLaplace(const int &n, const int &i, const int &j){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,i,j));
    sp.push_back(specialPoint(n,i-1,j));
    sp.push_back(specialPoint(n,i+1,j));
    sp.push_back(specialPoint(n,i,j-1));
    sp.push_back(specialPoint(n,i,j+1));
    sp.push_back(specialPoint(n,i+1,j-1));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    b(3) = -1;
    b(4) = -1;
    ColVector res1 = solve(A,b);

    sp[5] = specialPoint(n,i+1,j+1);
    A = specialGenA(sp);
    ColVector res2 = solve(A,b);

    sp[5] = specialPoint(n,i-1,j-1);
    A = specialGenA(sp);
    ColVector res3 = solve(A,b);

    sp[5] = specialPoint(n,i-1,j+1);
    A = specialGenA(sp);
    ColVector res4 = solve(A,b);

    ColVector res(9);
    for(int i = 0; i < 5; i++)
        res(i) = res1(i) + res2(i) + res3(i) + res4(i);
    res(5) = res1(5);
    res(6) = res2(5);
    res(7) = res3(5);
    res(8) = res4(5);
    return 0.25 * res;
}

ColVector specialSolveLeftNeumann(const int &n, const int &i, const int &j){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,i,j));
    sp.push_back(specialPoint(n,i+1,j));
    sp.push_back(specialPoint(n,i+2,j));
    sp.push_back(specialPoint(n,i,j-1));
    sp.push_back(specialPoint(n,i,j+1));
    sp.push_back(specialPoint(n,i+1,j+1));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    b(1) = -1;
    ColVector res1 = solve(A,b);

    sp[5] = specialPoint(n,i+1,j-1);
    A = specialGenA(sp);
    ColVector res2 = solve(A,b);

    ColVector res(7);
    for(int i = 0; i < 5; i++)
        res(i) = res1(i) + res2(i);
    res(5) = res1(5);
    res(6) = res2(5);
    return 0.5 * res;
}

ColVector specialSolveRightNeumann(const int &n, const int &i, const int &j){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,i,j));
    sp.push_back(specialPoint(n,i-1,j));
    sp.push_back(specialPoint(n,i-2,j));
    sp.push_back(specialPoint(n,i,j-1));
    sp.push_back(specialPoint(n,i,j+1));
    sp.push_back(specialPoint(n,i-1,j+1));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    b(1) = 1;
    ColVector res1 = solve(A,b);

    sp[5] = specialPoint(n,i-1,j-1);
    A = specialGenA(sp);
    ColVector res2 = solve(A,b);

    ColVector res(7);
    for(int i = 0; i < 5; i++)
        res(i) = res1(i) + res2(i);
    res(5) = res1(5);
    res(6) = res2(5);
    return 0.5 * res;
}

ColVector specialSolveUpNeumann(const int &n, const int &i, const int &j){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,i,j));
    sp.push_back(specialPoint(n,i,j-1));
    sp.push_back(specialPoint(n,i,j-2));
    Matrix A(3,3);
    for(int i = 0; i < 3; i++){
        A(0,i) = 1;
        A(1,i) = sp[i].y - sp[0].y;
        A(2,i) = 0.5 * (sp[i].y - sp[0].y) * (sp[i].y - sp[0].y);
    }

    ColVector b(3);
    b(1) = 1;
    return solve(A,b);
}

ColVector specialSolveLowNeumann(const int &n, const int &i, const int &j){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,i,j));
    sp.push_back(specialPoint(n,i,j+1));
    sp.push_back(specialPoint(n,i,j+2));
    sp.push_back(specialPoint(n,i-1,j));
    sp.push_back(specialPoint(n,i+1,j));
    sp.push_back(specialPoint(n,i-1,j+1));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    ColVector normV = normalVector(1.0*i/n);
    b(1) = normV(0);
    b(2) = normV(1);
    ColVector res1 = solve(A,b);

    sp[5] = specialPoint(n,i+1,j+1);
    A = specialGenA(sp);
    ColVector res2 = solve(A,b);

    ColVector res(7);
    for(int i = 0; i < 5; i++)
        res(i) = res1(i) + res2(i);
    res(5) = res1(5);
    res(6) = res2(5);
    return 0.5 * res;
}

ColVector specialSolveLowNeumann_LDcorner(const int &n){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,0,0));
    sp.push_back(specialPoint(n,0,1));
    sp.push_back(specialPoint(n,0,2));
    sp.push_back(specialPoint(n,1,0));
    sp.push_back(specialPoint(n,1,1));
    sp.push_back(specialPoint(n,2,0));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    ColVector normV = normalVector(0.0);
    b(1) = normV(0);
    b(2) = normV(1);
    return solve(A,b);
}

ColVector specialSolveLowNeumann_RDcorner(const int &n){
    std::vector<specialPoint> sp;
    sp.push_back(specialPoint(n,n,0));
    sp.push_back(specialPoint(n,n,1));
    sp.push_back(specialPoint(n,n,2));
    sp.push_back(specialPoint(n,n-1,0));
    sp.push_back(specialPoint(n,n-1,1));
    sp.push_back(specialPoint(n,n-2,0));
    Matrix A = specialGenA(sp);
    ColVector b(6);
    ColVector normV = normalVector(1.0);
    b(1) = normV(0);
    b(2) = normV(1);
    return solve(A,b);
}

//--------------------------------Solver 2D--------------------------------

Solver2D::Solver2D(const Operator & restriction, const Operator & prolongation)
    : rootSolver(restriction, prolongation){
    dim = 2;
    nineStencil = false;
}

void Solver2D::init(const int &_n, const Function &f, const Function &g, const std::string &bon, const BType &_bonDetail){
    n = _n;
    bondary = bon;
    bonDetail = _bonDetail;
    cur = zeroGridCol(n,2);
    b = zeroGridCol(n,2);
    double h = 1.0/n;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            if(irregular){
                specialPoint p(n,i,j);
                b(IDX(n,i,j)) = f(p.x, p.y);
            } else if(nineStencil){
                b(IDX(n,i,j)) = f(i*h, j*h) + h*h/12.0 * f.delta(i*h, j*h);
            } else {
                b(IDX(n,i,j)) = f(i*h, j*h);
            }
        }
    for(int i = 0; i <= n; i++){
        if(irregular){
            b(IDX(n,i,0)) = g(i*h, lowBon(i*h));
        } else {
            b(IDX(n,i,0)) = g(i*h, 0);
        }
        b(IDX(n,i,n)) = g(i*h, 1);
        if(i==0 || i==n) continue;
        b(IDX(n,0,i)) = g(0, i*h);
        b(IDX(n,n,i)) = g(1, i*h);
    }
    if(bon=="Neumann" || bon=="mixed" && bonDetail("x=0")=="Neumann" && bonDetail("x=1")=="Neumann"
       && bonDetail("Down Bondary")=="Neumann" && bonDetail("y=1")=="Neumann")
        pure_Neumann = true;
}

SparseMatrix Solver2D::getAh(const int &n) const{
    double h = 1.0/n;
    std::vector<Triple> ele;
    for(int i = 1; i < n; i++)
        for(int j = 1; j < n; j++){
            if(irregular){
                ColVector coef = specialSolveLaplace(n,i,j);
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j), coef(0)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j), coef(1)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j), coef(2)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j-1), coef(3)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j+1), coef(4)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j-1), coef(5)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j+1), coef(6)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j-1), coef(7)));
                ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j+1), coef(8)));
            } else {
                if(nineStencil){
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j), 10.0/(3.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j), -2.0/(3.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j), -2.0/(3.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j-1), -2.0/(3.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j+1), -2.0/(3.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j-1), -1.0/(6.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j-1), -1.0/(6.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j+1), -1.0/(6.0*h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j+1), -1.0/(6.0*h*h)));
                } else {
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j), 4.0/(h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i-1,j), -1.0/(h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i+1,j), -1.0/(h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j-1), -1.0/(h*h)));
                    ele.emplace_back(Triple(IDX(n,i,j), IDX(n,i,j+1), -1.0/(h*h)));
                }
            }
        }
    for(int i = 0; i <= n; i++){
        // bondary y=0 or irregular low-boundary
        if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("Down Bondary")=="Dirichlet"){
            ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,0), 1.0));
        } else {
            if(irregular){
                if(i==0){
                    ColVector coef = specialSolveLowNeumann_LDcorner(n);
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,0), coef(0)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,1), coef(1)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,2), coef(2)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i+1,0), coef(3)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i+1,1), coef(4)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i+2,0), coef(5)));
                } else if(i==n){
                    ColVector coef = specialSolveLowNeumann_RDcorner(n);
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n,0), coef(0)));
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n,1), coef(1)));
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n,2), coef(2)));
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n-1,0), coef(3)));
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n-1,1), coef(4)));
                    ele.emplace_back(Triple(IDX(n,n,0), IDX(n,n-2,0), coef(5)));
                } else {
                    ColVector coef = specialSolveLowNeumann(n, i, 0);
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,0), coef(0)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,1), coef(1)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,2), coef(2)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i-1,0), coef(3)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i+1,0), coef(4)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i-1,1), coef(5)));
                    ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i+1,1), coef(6)));
                }
            } else {
                ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,0), 1.5/h));
                ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,1), -2.0/h));
                ele.emplace_back(Triple(IDX(n,i,0), IDX(n,i,2), 0.5/h));
            }
        }
        // bondary y=1
        if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("y=1")=="Dirichlet"){
            ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n), 1.0));
        } else {
            if(irregular){
                ColVector coef = specialSolveUpNeumann(n,i,n);
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n), coef(0)));
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n-1), coef(1)));
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n-2), coef(2)));
            } else {
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n), 1.5/h));
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n-1), -2.0/h));
                ele.emplace_back(Triple(IDX(n,i,n), IDX(n,i,n-2), 0.5/h));
            }
        }
        if(i==0 || i==n) continue;
        // bondary x=0
        if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("x=0")=="Dirichlet"){
            ele.emplace_back(Triple(IDX(n,0,i), IDX(n,0,i), 1.0));
        } else {
            if(irregular){
                ColVector coef = specialSolveLeftNeumann(n, 0, i);
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,0,i), coef(0)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,1,i), coef(1)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,2,i), coef(2)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,0,i-1), coef(3)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,0,i+1), coef(4)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,1,i+1), coef(5)));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,1,i-1), coef(6)));
            } else {
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,0,i), 1.5/h));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,1,i), -2.0/h));
                ele.emplace_back(Triple(IDX(n,0,i), IDX(n,2,i), 0.5/h));
            }
        }
        // bondary x=1
        if(bondary == "Dirichlet" || bondary == "mixed" && bonDetail("x=1")=="Dirichlet"){
            ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n,i), 1.0));
        } else {
            if(irregular){
                ColVector coef = specialSolveRightNeumann(n, n, i);
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n,i), coef(0)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-1,i), coef(1)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-2,i), coef(2)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n,i-1), coef(3)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n,i+1), coef(4)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-1,i+1), coef(5)));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-1,i-1), coef(6)));
            } else {
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n,i), 1.5/h));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-1,i), -2.0/h));
                ele.emplace_back(Triple(IDX(n,n,i), IDX(n,n-2,i), 0.5/h));
            }
        }
    }
    return SparseMatrix( (n+1)*(n+1), (n+1)*(n+1), ele );
}

void Solver2D::output(std::ostream &out) const{
    out << std::fixed << std::setprecision(16);
    double h = 1.0/n;
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++){
            double x = irregular ? specialPoint(n,i,j).x : i*h;
            double y = irregular ? specialPoint(n,i,j).y : j*h;
            out << x << " " << y << " " << cur(IDX(n,i,j)) << "\n";
        }
}

double Solver2D::calcNeumannC(const Function &u){
    double l = 1e100, r = -1e100;
    double h = 1.0/n;
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++){
            double x = irregular ? specialPoint(n,i,j).x : i*h;
            double y = irregular ? specialPoint(n,i,j).y : j*h;
            l = std::min(l, u(x, y) - cur(IDX(n,i,j)));
            r = std::max(r, u(x, y) - cur(IDX(n,i,j)));
        }
    return neumannC = (l+r)/2;
}

double Solver2D::checkError(const Function &u, const Norm &norm) const{
    std::vector<double> vals;
    double h = 1.0/n;
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= n; j++){
            double x = irregular ? specialPoint(n,i,j).x : i*h;
            double y = irregular ? specialPoint(n,i,j).y : j*h;
            vals.emplace_back( cur(IDX(n,i,j)) + neumannC - u(x, y) );
        }
    return norm(vals);
}