#include "solver.h"
#include "matrix.h"
#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
#include <algorithm>

const double _epsL = 1e-12;

//-----------------------------Circle----------------------------------

Circle::Circle()
    : Circle(0, 0, 0){}

Circle::Circle(const double &_x, const double &_y, const double &_R)
    : cx(_x), cy(_y), R(_R){}

bool Circle::inCircle(const double &x, const double &y) const{
    return (x-cx)*(x-cx) + (y-cy)*(y-cy) < R*R + _epsL;
}

bool Circle::onCircle(const double &x, const double &y) const{
    return fabs( (x-cx)*(x-cx) + (y-cy)*(y-cy) - R*R ) < _epsL;
}

double Circle::crossX(const double &x0, const int &i) const{
    double dy;
    if(fabs(x0-cx-R) < _epsL || fabs(x0-cx+R) < _epsL) dy = 0;
    else dy = sqrt(R*R - (x0-cx)*(x0-cx));
    return (i==0) ? cy-dy : cy+dy;
}

double Circle::crossY(const double &y0, const int &i) const{
    double dx;
    if(fabs(y0-cy-R) < _epsL || fabs(y0-cy+R) < _epsL) dx = 0;
    else dx = sqrt(R*R - (y0-cy)*(y0-cy));
    return (i==0) ? cx-dx : cx+dx;
}

//------------------------------Point------------------------------------

int Point::getID() const{
    return pid;
}

double Point::getX() const{
    return x;
}

double Point::getY() const{
    return y;
}

//-----------------------------Function2D----------------------------------

ColVector Function2D::grad(const double &x, const double &y) const{
    ColVector res(2);
    res(1) = ((*this)(x+_epsL,y) - (*this)(x-_epsL, y)) / (2*_epsL);
    res(2) = ((*this)(x,y+_epsL) - (*this)(x, y-_epsL)) / (2*_epsL);
    return res;
}

double Function2D::div(const double &x, const double &y) const{
    return sum(grad(x, y));
}

double Function2D::delta(const double &x, const double &y) const{
    return ( (*this)(x+_epsL,y) + (*this)(x-_epsL,y) + (*this)(x,y+_epsL) + (*this)(x,y-_epsL) - 4*(*this)(x,y) ) / (_epsL * _epsL);
}

//-----------------------------BType----------------------------------

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

//-----------------------------Solver----------------------------------

int Solver::P(const int &x, const int &y) const{
    return x * (m+2) + y;
}

int Solver::irpID(const double &x, const double &y) const{
    for(const Point &p : irp)
        if( fabs(p.getX()-x)<_epsL && fabs(p.getY()-y)<_epsL )
            return p.getID();
    std::cerr << "[Error] There's no irregular point at (" << x << "," << y << ")" << std::endl;
    return -1;
}

bool Solver::isIrp(const double &x, const double &y) const{
    for(const Point &p : irp)
        if( fabs(p.getX()-x)<_epsL && fabs(p.getY()-y)<_epsL )
            return true;
    return false;
}

Point Solver::toPoint(const int &i, const int &j) const{
    return Point(P(i,j), 1.0*i/(m+1), 1.0*j/(m+1));
}

void Solver::solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, 
                const double &cx, const double &cy, const double &R){
    solve(f, g, _m, bondary, cx, cy, R, BType());
};

void Solver::solve(Function2D &f, Function2D &g, const int &_m, std::string bondary){
    solve(f, g, _m, bondary, BType());
}

void Solver::solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, const BType &bondaryDetail){
    solve(f, g, _m, bondary, 0, 0, 0, bondaryDetail);
}

double sqr(const double &x){
    return x * x;
}

// 添加一个 irregular point 的 Neumann 边值条件下的方程
void irNeumannEquation(Matrix &A, const int &cnt, const double &cosa, const double &sina, const std::vector<Point> &pnt){
    Matrix coef(6);
    ColVector b(6);
    for(int i = 0; i < 6; i++){
        double dx = pnt[i].getX() - pnt[0].getX();
        double dy = pnt[i].getY() - pnt[0].getY();
        coef(0,i) = 1;
        coef(1,i) = dx;
        coef(2,i) = dy;
        coef(3,i) = 0.5*dx*dx;
        coef(4,i) = 0.5*dy*dy;
        coef(5,i) = dx*dy;
    }
    b(1) = cosa;
    b(2) = sina;
    ColVector res = coef.solveDestructiveness(b);
    for(int i = 0; i < 6; i++)
        A(cnt, pnt[i].getID()) = res(i);
}

void Solver::solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, 
                const double &cx, const double &cy, const double &R, const BType &bondaryDetail){
    m = _m;
    double h = 1.0/(m+1), x, y, lx, ly, rx, ry, a;
    int pcnt = (m+2)*(m+2);
    c = Circle(cx, cy, R);
    int firstBdryCond;

    // Find all irregular points
    for(int i = 1; i <= m; i++){
        x = i * h;
        if(fabs(x-cx) <= R + _epsL){
            Point a = Point(pcnt, x, c.crossX(x,0));
            Point b = Point(pcnt+1, x, c.crossX(x,1));
            if(fabs(a.getY()-b.getY()) < _epsL){
                if(fabs(a.getY() - h*floor(a.getY()/h)) < _epsL)
                    irp.push_back(a), pcnt++;
            }
            else{
                irp.push_back(a);
                irp.push_back(b);
                pcnt += 2;
            }
        }
    }
    for(int j = 1; j <= m; j++){
        y = j * h;
        if(fabs(y-cy) <= R + _epsL){
            Point a = Point(pcnt, c.crossY(y,0), y);
            Point b = Point(pcnt+1, c.crossY(y,1), y);
            if(fabs(a.getX()-b.getX()) < _epsL){
                if(fabs(a.getX() - h*floor(a.getX()/h)) < _epsL){
                    irp.push_back(a);
                    pcnt++;
                }
            } else {
                if(!isIrp(a.getX(),a.getY())){
                    irp.push_back(a);
                    pcnt++;
                }
                if(!isIrp(b.getX(),b.getY())){
                    irp.push_back(b);
                    pcnt++;
                }
            }
        }
    }
    Matrix A(pcnt, pcnt);
    ColVector b(pcnt);

    // Equation
    int cnt = 0;
    for(int i = 1, irpid; i <= m; i++)
        for(int j = 1; j <= m; j++){
            x = 1.0*i/(m+1);
            y = 1.0*j/(m+1);
            if(c.inCircle(x,y)){
                A(cnt, P(i,j)) = 1;
                b(cnt++) = 0;
                continue;
            }
            lx =  x - h;
            rx =  x + h;
            ly =  y - h;
            ry =  y + h;
            b(cnt) = h*h * f( 1.0*i/(m+1), 1.0*j/(m+1) );
            if(!c.inCircle(lx,y) && !c.inCircle(rx,y)){
                A(cnt, P(i+1,j)) = -1;
                A(cnt, P(i-1,j)) = -1;
                A(cnt, P(i,j)) = 2;
            } else if(c.inCircle(lx,y)){
                a = (x - c.crossY(y,1)) / h;
                irpid = irpID(x-a*h, y);
                A(cnt, P(i,j)) = 2/a;
                A(cnt, P(i+1,j)) = -2/(1+a);
                A(cnt, irpid) = -2/(a*(1+a));
            } else {
                a = (c.crossY(y,0) - x) / h;
                irpid = irpID(x+a*h, y);
                A(cnt, P(i,j)) = 2/a;
                A(cnt, P(i-1,j)) = -2/(1+a);
                A(cnt, irpid) = -2/(a*(1+a));
            }

            if(!c.inCircle(x,ly) && !c.inCircle(x,ry)){
                A(cnt, P(i,j+1)) = -1;
                A(cnt, P(i,j-1)) = -1;
                A(cnt, P(i,j)) += 2;
            } else if(c.inCircle(x,ly)){
                a = (y - c.crossX(x,1)) / h;
                irpid = irpID(x, y-a*h);
                A(cnt, P(i,j)) += 2/a;
                A(cnt, P(i,j+1)) = -2/(1+a);
                A(cnt, irpid) = -2/(a*(1+a));
            } else {
                a = (c.crossX(x,0) - y) / h;
                irpid = irpID(x, y+a*h);
                A(cnt, P(i,j)) += 2/a;
                A(cnt, P(i,j-1)) = -2/(1+a);
                A(cnt, irpid) = -2/(a*(1+a));
            }
            cnt++;
        }

    if(bondary == "Dirichlet"){
        // Dirichlet Bondary Conditions
        pure_Neumann = false;
        for(int i = 0; i <= m+1; i++){
            A(cnt, P(i,0)) = 1;
            b(cnt++) = g(1.0*i/(m+1), 0);
            A(cnt, P(i,m+1)) = 1;
            b(cnt++) = g(1.0*i/(m+1), 1);
            if(i==0 || i==m+1) continue;
            A(cnt, P(0,i)) = 1;
            b(cnt++) = g(0, 1.0*i/(m+1));
            A(cnt, P(m+1,i)) = 1;
            b(cnt++) = g(1, 1.0*i/(m+1));
        }
        for(const Point &p : irp){
            A(cnt, p.getID()) = 1;
            b(cnt++) = g(p.getX(), p.getY());
        }
    } else if(bondary == "Neumann"){
        // Neumann Bondary Conditions
        pure_Neumann = true;
        A(cnt, P(0,0))++;
        std::cout << "[Warning] The solution with Neumann condition will have a CONSTANT difference from the real solution." << std::endl;
        for(int i = 0; i <= m+1; i++){
            A(cnt, P(i,0)) = 1.5;
            A(cnt, P(i,1)) = -2;
            A(cnt, P(i,2)) = 0.5;
            b(cnt++) = g(1.0*i/(m+1), 0) * h;
            A(cnt, P(i,m+1)) = 1.5;
            A(cnt, P(i,m)) = -2;
            A(cnt, P(i,m-1)) = 0.5;
            b(cnt++) = g(1.0*i/(m+1), 1) * h;
            if(i==0 || i==m+1) continue;
            A(cnt, P(0,i)) = 1.5;
            A(cnt, P(1,i)) = -2;
            A(cnt, P(2,i)) = 0.5;
            b(cnt++) = g(0, 1.0*i/(m+1)) * h;
            A(cnt, P(m+1,i)) = 1.5;
            A(cnt, P(m,i)) = -2;
            A(cnt, P(m-1,i)) = 0.5;
            b(cnt++) = g(1, 1.0*i/(m+1)) * h;
        }
    } else if(bondary == "mixed"){
        firstBdryCond = cnt;
        pure_Neumann = true;
        if(bondaryDetail("x=0")=="Dirichlet"){
            for(int i = 0; i <= m+1; i++){
                A(cnt, P(0,i)) = 1;
                b(cnt++) = g(0, 1.0*i/(m+1));
            }
            pure_Neumann = false;
        } else {
            for(int i = 0; i <= m+1; i++){
                if(c.inCircle(0,i*h)){
                    A(cnt, P(0,i)) = 1;
                    b(cnt++) = 0;
                    continue;
                }
                A(cnt, P(0,i)) = 1.5;
                A(cnt, P(1,i)) = -2;
                A(cnt, P(2,i)) = 0.5;
                b(cnt++) = g(0, 1.0*i/(m+1)) * h;
            }
        }

        if(bondaryDetail("x=1")=="Dirichlet"){
            for(int i = 0; i <= m+1; i++){
                A(cnt, P(m+1,i)) = 1;
                b(cnt++) = g(1, 1.0*i/(m+1));
            }
            pure_Neumann = false;
        } else {
            for(int i = 0; i <= m+1; i++){
                if(c.inCircle(1,i*h)){
                    A(cnt, P(m+1,i)) = 1;
                    b(cnt++) = 0;
                    continue;
                }
                A(cnt, P(m+1,i)) = 1.5;
                A(cnt, P(m,i)) = -2;
                A(cnt, P(m-1,i)) = 0.5;
                b(cnt++) = g(1, 1.0*i/(m+1)) * h;
            }
        }

        if(bondaryDetail("y=0")=="Dirichlet"){
            for(int i = 1; i <= m; i++){
                A(cnt, P(i,0)) = 1;
                b(cnt++) = g(1.0*i/(m+1),0);
            }
            pure_Neumann = false;
        } else {
            for(int i = 1; i <= m; i++){
                if(c.inCircle(i*h,0)){
                    A(cnt, P(i,0)) = 1;
                    b(cnt++) = 0;
                    continue;
                }
                A(cnt, P(i,0)) = 1.5;
                A(cnt, P(i,1)) = -2;
                A(cnt, P(i,2)) = 0.5;
                b(cnt++) = g(1.0*i/(m+1),0) * h;
            }
        }

        if(bondaryDetail("y=1")=="Dirichlet"){
            for(int i = 1; i <= m; i++){
                A(cnt, P(i,m+1)) = 1;
                b(cnt++) = g(1.0*i/(m+1),1);
            }
            pure_Neumann = false;
        } else {
            for(int i = 1; i <= m; i++){
                if(c.inCircle(i*h,1)){
                    A(cnt, P(i,m+1)) = 1;
                    b(cnt++) = 0;
                    continue;
                }
                A(cnt, P(i,m+1)) = 1.5;
                A(cnt, P(i,m)) = -2;
                A(cnt, P(i,m-1)) = 0.5;
                b(cnt++) = g(1.0*i/(m+1),1) * h;
            }
        }

        if(R && bondaryDetail("D")=="Dirichlet"){
            for(const Point &p : irp){
                A(cnt, p.getID()) = 1;
                b(cnt++) = g(p.getX(), p.getY());
            }
            pure_Neumann = false;
        }
    } else{
        std::cerr << "[Error] There's no bondary condition named '" << bondary << "', try 'Dirichlet', 'Neumann' or 'mixeded'." << std::endl;
        exit(-1);
    }

    if(bondary=="Neumann" || bondary=="mixed" && R && bondaryDetail("D")=="Neumann"){
        for(const Point &p : irp){
            std::vector<Point> pnt(6);
            double cosb = (p.getX()-cx)/R;
            double sinb = (p.getY()-cy)/R;
            x = p.getX();
            y = p.getY();
            pnt[0] = p;
            int i = (int)round(x/h), j = (int)round(y/h);
            if( fabs(x-h*i)<_epsL && fabs(y-h*j)<_epsL ){ //在网格点上
                int st1 = (x-cx>=0) ? 1 : -1;
                int st2 = (y-cy>=0) ? 1 : -1;
                pnt[1] = toPoint(i+st1, j);
                pnt[2] = toPoint(i+2*st1, j);
                pnt[3] = toPoint(i+st1, j+st2);
                pnt[4] = toPoint(i, j+st2);
                pnt[5] = toPoint(i, j+2*st2);
            } else if( fabs(x-h*i)<_epsL && x-cx>=0 && y-cy>=0 ){ //在网格某条竖线上，在圆右上角
                int i = (int)round(x/h), j = (int)ceil(y/h);
                pnt[1] = toPoint(i,j+1);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i+1,j+1);
                pnt[4] = toPoint(i+1,j);
                pnt[5] = toPoint(i+2,j);
            } else if( fabs(y-h*j)<_epsL && x-cx>=0 && y-cy>=0 ){ //在网格某条横线上，在圆右上角
                int i = (int)ceil(x/h), j = (int)round(y/h);
                pnt[1] = toPoint(i+1,j);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i+1,j+1);
                pnt[4] = toPoint(i,j+1);
                pnt[5] = toPoint(i,j+2);
            } else if( fabs(x-h*i)<_epsL && x-cx>=0 && y-cy<0 ){ //在网格某条竖线上，在圆右下角
                int i = (int)round(x/h), j = (int)floor(y/h);
                pnt[1] = toPoint(i,j-1);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i+1,j-1);
                pnt[4] = toPoint(i+1,j);
                pnt[5] = toPoint(i+2,j);
            } else if( fabs(y-h*j)<_epsL && x-cx>=0 && y-cy<0 ){ //在网格某条横线上，在圆右下角
                int i = (int)ceil(x/h), j = (int)round(y/h);
                pnt[1] = toPoint(i+1,j);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i+1,j-1);
                pnt[4] = toPoint(i,j-1);
                pnt[5] = toPoint(i,j-2);
            } else if( fabs(x-h*i)<_epsL && x-cx<0 && y-cy>=0 ){ //在网格某条竖线上，在圆左上角
                int i = (int)round(x/h), j = (int)ceil(y/h);
                pnt[1] = toPoint(i,j+1);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i-1,j+1);
                pnt[4] = toPoint(i-1,j);
                pnt[5] = toPoint(i-2,j);
            } else if( fabs(y-h*j)<_epsL && x-cx<0 && y-cy>=0 ){ //在网格某条横线上，在圆左上角
                int i = (int)floor(x/h), j = (int)round(y/h);
                pnt[1] = toPoint(i-1,j);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i-1,j+1);
                pnt[4] = toPoint(i,j+1);
                pnt[5] = toPoint(i,j+2);
            } else if( fabs(x-h*i)<_epsL && x-cx<0 && y-cy<0 ){ //在网格某条竖线上，在圆左下角
                int i = (int)round(x/h), j = (int)floor(y/h);
                pnt[1] = toPoint(i,j-1);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i-1,j-1);
                pnt[4] = toPoint(i-1,j);
                pnt[5] = toPoint(i-2,j);
            } else if( fabs(y-h*j)<_epsL && x-cx<0 && y-cy<0 ){ //在网格某条横线上，在圆左下角
                int i = (int)floor(x/h), j = (int)round(y/h);
                pnt[1] = toPoint(i-1,j);
                pnt[2] = toPoint(i,j);
                pnt[3] = toPoint(i-1,j-1);
                pnt[4] = toPoint(i,j-1);
                pnt[5] = toPoint(i,j-2);
            }
            irNeumannEquation(A, cnt, cosb, sinb, pnt);
            b(cnt++) = g(x,y);
        }
    }

    if(bondary == "mixed" && pure_Neumann){
        A(firstBdryCond, P(0,0))++;
        std::cout << "[Warning] Your mixed bondary is actually pure-Neumann." << std::endl;
        std::cout << "[Warning] The solution with Neumann condition will have a CONSTANT difference from the real solution." << std::endl;
    }

    // Solve the linear system
    Uval = A.solveDestructiveness(b);
}

bool Solver::inRange(const double &x, const double &y) const{
    return x>0 && x<1 && y>0 && y<1 && !c.inCircle(x,y);
}

double Solver::operator () (const double &x, const double &y) const{
    if(m==0){
        std::cerr << "[Error] Cannot call operator () in an empty Solver!" << std::endl;
        exit(-1);
    } else if(!inRange(x,y)){
        std::cerr << "[Error] Operator () :: (" << x << "," << y << ") out of range!" << std::endl;
        exit(-1);
    }
    int i = (int) float(x * (m + 1));
    int j = (int) float(y * (m + 1));
    double a = x*(m+1) - i;
    double b = y*(m+1) - j;
    return (1-a)*(1-b)*Uval(P(i,j)) + (1-a)*b*Uval(P(i,j+1)) + a*(1-b)*Uval(P(i+1,j)) + a*b*Uval(P(i+1,j+1));
}

std::ostream & operator << (std::ostream & out, const Solver & ps){
    typedef std::pair< std::pair<double,double>, double > Vp;
    std::vector<Vp> vec;
    out << std::fixed << std::setprecision(12); 
    for(int i = 1; i <= ps.m; i++)
        for(int j = 1; j <= ps.m; j++){
            double x = 1.0*i/(ps.m+1);
            double y = 1.0*j/(ps.m+1);
            if(!ps.inRange(x,y)) continue;
            vec.push_back( std::make_pair( std::make_pair(x,y), ps.Uval(ps.P(i,j)) ) );
        }
    for(const Point & p : ps.irp)
        vec.push_back( std::make_pair( std::make_pair(p.getX(),p.getY()), ps.Uval(p.getID())) );
    std::sort(vec.begin(), vec.end());
    for(const Vp & p : vec)
        out << p.first.first << " " << p.first.second << " " << p.second << "\n";
    return out;
}

bool Solver::isPureNeumann() const{
    return pure_Neumann;
}

void Solver::calcNeumannC(Function2D &u){
    have_Neumann_C = true;
    double l = 1e100, r = -1e100;
    for(int i = 1; i <= m; i++)
        for(int j = 1; j <= m; j++){
            double x = 1.0 * i / (m+1);
            double y = 1.0 * j / (m+1);
            if(!inRange(x,y)) continue;
            l = std::min( l, u(x,y) - Uval(P(i,j)) );
            r = std::max( r, u(x,y) - Uval(P(i,j)) );
        }
    for(const Point &p : irp){
        r = std::max( r, u(p.getX(),p.getY()) - Uval(p.getID()) );
        l = std::min( l, u(p.getX(),p.getY()) - Uval(p.getID()) );
    }
    Neumann_C = (l+r)/2;
    std::cerr << "The numerical solutions should add a constant " << Neumann_C << std::endl;
}

double Solver::checkError(Function2D &u, const Norm &norm){
    if(m==0){
        std::cerr << "[Error] Cannot check the error of an empty solver" << std::endl;
        exit(-1);
    }
    if(pure_Neumann && !have_Neumann_C){
        calcNeumannC(u);
    }
    std::vector<double> vec;
    for(int i = 1; i <= m; i++)
        for(int j = 1; j <= m; j++){
            double x = 1.0 * i / (m+1);
            double y = 1.0 * j / (m+1);
            if(!inRange(x,y)) continue;
            vec.push_back(Uval(P(i,j)) + Neumann_C - u(x,y));
        }
    for(const Point &p : irp){
        vec.push_back( Uval(p.getID()) + Neumann_C - u(p.getX(),p.getY()) );
    }
    return norm(vec);
}
