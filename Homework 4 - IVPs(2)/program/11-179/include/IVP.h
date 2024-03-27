#ifndef _IVP_H_
#define _IVP_H_

#include "matrix.h"
#include <cstring>
#include <vector>

class TimeFunction{
public:
    virtual ColVector operator () (const ColVector &x, const double &t) const = 0;
};


class TimeIntegrator{
protected:
    std::vector<ColVector> dsol;
    std::vector<ColVector> sol;
    double maxTime;

public:
    TimeIntegrator(): maxTime(-1){};
    virtual void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize) = 0;
    virtual void solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps) = 0;
    virtual ColVector at(const double &t) const = 0;    //返回t时刻的求解结果，如果t不是求解的时间结点，则以三阶Hermite插值给出
    virtual void output(const std::string &fname) const = 0;       //输出到文件，每行格式：t x1 x2 ... xn
    virtual void denseOutput(const std::string &fname) const = 0;
    void denseDiscreteOutput(const std::string &fname, const double &step) const;
};


class TimeIntegrator_ConstStep : public TimeIntegrator{  //常步长求解器基类
protected:
    double timeStep;

public:
    ColVector at(const double &t) const;
    void output(const std::string &fname) const;
    void denseOutput(const std::string &fname) const;
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps);
};


class TimeIntegrator_VariativeStep : public TimeIntegrator{  //变步长求解器基类
protected:
    std::vector<double> timePoint;

public:
    ColVector at(const double &t) const;
    void output(const std::string &fname) const;
    void denseOutput(const std::string &fname) const;
    void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
};

#endif