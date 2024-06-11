#ifndef _IVP_H_
#define _IVP_H_

#include "matrix.h"
#include <cstring>
#include <vector>
#include <map>

class TimeFunction{
protected:
    bool __isLinear;
public:
    TimeFunction(): __isLinear(false){}
    ColVector operator () (const ColVector &x, const double &t) const{
        ColVector res = x;
        compute(x, t, res);
        return res;
    }
    virtual ColVector dp(const ColVector &p, const ColVector &q) const {return ColVector();}
    virtual ColVector dq(const ColVector &p, const ColVector &q) const {return ColVector();}
    virtual ColVector solve(const Matrix &A, const ColVector &c, const ColVector &U0, const double &t0, const double &k) const {return ColVector();}
    bool isLinear() const{return __isLinear;}
    virtual void compute(const ColVector &x, const double &t, ColVector &res) const = 0;
};


class TimeIntegrator{
protected:
    std::vector<ColVector> sol, dsol;
    double maxTime, maxStep;
    int maxorder;

public:
    TimeIntegrator(): maxTime(-1), maxStep(-1){};
    void solveWithInfo(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize);
    void solveWithInfo(TimeFunction &f, const ColVector &x0, const double &T, const double &eps);
    void setMaxStep(const double &t);
    virtual void solve(TimeFunction &f, const ColVector &x0, const double &T, const int &gridSize) = 0;
    virtual void solve(TimeFunction &f, const ColVector &x0, const double &T, const double &eps) = 0;
    ColVector getSolByID(const int &i) const;                                        //返回第i个时刻的解
    virtual ColVector at(const double &t) const = 0;                                 //返回t时刻的求解结果，如果t不是求解的时间结点，则以三阶Hermite插值给出
    virtual void output(const std::string &fname) const = 0;                         //输出数值解到文件，每行格式：t x1 x2 ... xn
    virtual void denseOutput(const std::string &fname) const = 0;                    //输出对数值解做三次二阶Hermite插值的结果
    void denseDiscreteOutput(const std::string &fname, const double &step) const;    //借助denseOutput的结果每隔step步长输出一个点值
};


class TimeIntegratorFactory{
public:
    using CreateTimeIntegratorCallBack = TimeIntegrator* (*)(int);
private:
    using CallbackMap = std::map<std::string, CreateTimeIntegratorCallBack>;
public:
    bool registerTimeIntegrator(const std::string &ID, CreateTimeIntegratorCallBack createFn);
    bool unregisterTimeIntegrator(const std::string &ID);
    TimeIntegrator* createTimeIntegrator(const std::string &ID, const int &order);
private:
    CallbackMap callbacks_;
private:
    TimeIntegratorFactory() = default;
    TimeIntegratorFactory(const TimeIntegratorFactory&) = default;
    TimeIntegratorFactory& operator = (const TimeIntegratorFactory&) = default;
    ~TimeIntegratorFactory() = default;
public:
    static TimeIntegratorFactory& Instance();
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
