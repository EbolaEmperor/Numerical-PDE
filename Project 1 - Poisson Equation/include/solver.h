#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "matrix.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>

class Circle{
private:
    double cx, cy, R;
public:
    Circle();
    Circle(const double &_x, const double &_y, const double &_R);
    bool inCircle(const double &x, const double &y) const;
    bool onCircle(const double &x, const double &y) const;
    double crossX(const double &x0, const int &i) const;
    double crossY(const double &y0, const int &i) const;
};

class Function2D{
public:
    virtual double operator () (const double &x, const double &y) const = 0;
    ColVector grad(const double &x, const double &y) const;
    double div(const double &x, const double &y) const;
    double delta(const double &x, const double &y) const;
};

class BType{
private:
    std::map<std::string, std::string> types;
public:
    void setBondary(const std::string& bon, const std::string& typ);
    std::string operator () (const std::string& bon) const;
};

class Point{
private:
    int pid;
    double x, y;
public:
    Point() {}
    Point(const int &pid, const double &x, const double &y):
        pid(pid), x(x), y(y) {}
    int getID() const;
    double getX() const;
    double getY() const;
};

class Norm{
public:
    virtual double operator () (const std::vector<double> &vec) const = 0;
};

class Solver{
private:
    int m;
    ColVector Uval;
    Circle c;
    std::vector<Point> irp;
    double Neumann_C;
    bool have_Neumann_C;
    bool pure_Neumann;

    int P(const int &x, const int &y) const;
    int irpID(const double &x, const double &y) const;
    bool isIrp(const double &x, const double &y) const;
    Point toPoint(const int &i, const int &j) const;

public:
    Solver(): m(0), Uval(ColVector()), have_Neumann_C(false), Neumann_C(0.0) {}
    void solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, 
                  const double &cx, const double &cy, const double &R);
    void solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, 
                  const double &cx, const double &cy, const double &R, const BType &bondaryDetail);
    void solve(Function2D &f, Function2D &g, const int &_m, std::string bondary);
    void solve(Function2D &f, Function2D &g, const int &_m, std::string bondary, const BType &bondaryDetail);

    bool inRange(const double &x, const double &y) const;
    double operator () (const double &x, const double &y) const;
    friend std::ostream & operator << (std::ostream & out, const Solver & ps);

    bool isPureNeumann() const;

    // 通过各网格点数值解与真解之间差值的平均，来确定纯Neumann条件下的常数C
    void calcNeumannC(Function2D &u);
    double checkError(Function2D &u, const Norm &norm);
};

#endif