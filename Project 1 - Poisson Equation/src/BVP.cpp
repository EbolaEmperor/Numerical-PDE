#include "BVP.h"
#include "solver.h"
#include "mathExpr.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <jsoncpp/json/json.h>
#include <string>
#include <vector>
#include <ctime>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
const double eps = 1e-15;

class Fun2D : public Function2D{
private:
    Expression e;
public:
    Fun2D(){
        e.addVariable("x");
        e.addVariable("y");
    }
    void setExpr(const std::string &expr){
        e.setExpr(expr);
    }
    double operator () (const double &x, const double &y) const{
        vector<double> parm;
        parm.push_back(x);
        parm.push_back(y);
        return e(parm);
    }
};

class BFun2D : public Function2D{
private:
    Expression e[5];
    Circle D;
public:
    BFun2D(){
        for(int i = 0; i < 5; i++){
            e[i].addVariable("x");
            e[i].addVariable("y");
        }
    }
    void setCircle(const double &x, const double &y, const double &R){
        D = Circle(x, y, R);
    }
    void addConstant(const std::string &bon, const std::string & name, const double &val){
        if(bon == "x=0") e[0].addConstant(name, val);
        else if(bon == "x=1") e[1].addConstant(name, val);
        else if(bon == "y=0") e[2].addConstant(name, val);
        else if(bon == "y=1") e[3].addConstant(name, val);
        else if(bon == "D") e[4].addConstant(name, val);
        else {
            cerr << "[Error] Unrecognized bondary '" << bon << "'" << endl;
            exit(-1);
        }
    }
    void setExpr(const std::string &bon, const std::string &expr){
        if(bon == "x=0") e[0].setExpr(expr);
        else if(bon == "x=1") e[1].setExpr(expr);
        else if(bon == "y=0") e[2].setExpr(expr);
        else if(bon == "y=1") e[3].setExpr(expr);
        else if(bon == "D") e[4].setExpr(expr);
        else {
            cerr << "[Error] Unrecognized bondary '" << bon << "'" << endl;
            exit(-1);
        }
    }
    double operator () (const double &x, const double &y) const{
        vector<double> parm;
        parm.push_back(x);
        parm.push_back(y);
        if(fabs(x) < eps) return e[0](parm);
        else if(fabs(x-1) < eps) return e[1](parm);
        else if(fabs(y) < eps) return e[2](parm);
        else if(fabs(y-1) < eps) return e[3](parm);
        else if(D.onCircle(x,y)) return e[4](parm);
        else{
            cerr << "[Error] Point (" << x << "," << y << ") is not on the bondary" << endl;
            exit(-1);
        } 
    }
};

class Norm_p : public Norm{
private:
    double p;
public:
    Norm_p(const int &p): p(p) {}
    double operator () (const std::vector<double> &vec) const{
        double sum = 0;
        for(const double &x : vec)
            sum += pow(fabs(x), p);
        sum /= vec.size();
        return pow(sum, 1.0/p);
    }
};

class Norm_inf : public Norm{
public:
    double operator () (const std::vector<double> &vec) const{
        double mx = 0;
        for(const double &x : vec)
            mx = std::max(mx, fabs(x));
        return mx;
    }
};

BVP::BVP(){}
BVP::BVP(const Json::Value & prob){
    read(prob);
}
BVP::BVP(const char * ifile){
    read(ifile);
}
BVP::BVP(const std::string & ifile){
    read(ifile);
}

void BVP::printProblem(){
    cout << "----------------------------------------------------------" << endl;
    cout << "BVP with " << problem["Condition Type"].asString() << " Condition" << endl;
    cout << "Reigeon Type: " << problem["Reigeon Type"].asString() << endl;
    if(problem["Reigeon Type"].asString() == "Irregular"){
        cout << "Irregular Circle: C=(" << problem["Center"][0].asDouble() << "," << problem["Center"][1].asDouble() << "), R=" << problem["Radius"].asDouble() << endl;
    }
    cout << "f(x,y) = " << problem["f"].asString() << endl;
    if( !problem["g"].isObject() ){
        cout << "g(x,y) = " << problem["g"].asString() << endl;
    } else {
        if(problem["Condition Type"].asString() == "mixed"){
            cout << "g(x,y) = " << problem["g"]["x=0"][0].asString() << "  if x=0  (" << problem["g"]["x=0"][1].asString() << ")" << endl;
            cout << "         " << problem["g"]["x=1"][0].asString() << "  if x=1  (" << problem["g"]["x=1"][1].asString() << ")" << endl;
            cout << "         " << problem["g"]["y=0"][0].asString() << "  if y=0  (" << problem["g"]["y=0"][1].asString() << ")" << endl;
            cout << "         " << problem["g"]["y=1"][0].asString() << "  if y=1  (" << problem["g"]["y=1"][1].asString() << ")" << endl;
            if(problem["Reigeon Type"].asString()=="Irregular")
                cout << "         " << problem["g"]["D"][0].asString() << "  if (x,y) in D  (" << problem["g"]["D"][1].asString() << ")" << endl;
        } else {
            cout << "g(x,y) = " << problem["g"]["x=0"].asString() << "  if x=0" << endl;
            cout << "         " << problem["g"]["x=1"].asString() << "  if x=1" << endl;
            cout << "         " << problem["g"]["y=0"].asString() << "  if y=0" << endl;
            cout << "         " << problem["g"]["y=1"].asString() << "  if y=1" << endl;
            if(problem["Reigeon Type"].asString()=="Irregular")
                cout << "         " << problem["g"]["D"].asString() << "  if (x,y) in D" << endl;
        }
    }
    cout << "Grid Size: " << problem["Grid Size"].asInt() << endl;
}

void BVP::read(const Json::Value &prob){
    problem = prob;
}

void BVP::read(const string & fname){
    Json::Reader reader;
    ifstream ifs;
    ifs.open(fname);
    if(!ifs.is_open()){
        cerr << "cannot read file " << fname << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }
}

void BVP::read(const char * fname){
    Json::Reader reader;
    ifstream ifs;
    ifs.open(fname);
    if(!ifs.is_open()){
        cerr << "cannot read file " << fname << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }
}

void BVP::solve(){
    Fun2D f, g1;
    BFun2D g2;
    Function2D *g;
    BType bonDtil;

    const int m = problem["Grid Size"].asInt();
    double cx = 0, cy = 0, R = 0;
    if(problem["Reigeon Type"].asString() == "Irregular"){
        cx = problem["Center"][0].asDouble();
        cy = problem["Center"][1].asDouble();
        R = problem["Radius"].asDouble();
    }

    f.setExpr(problem["f"].asString());
    if( !problem["g"].isObject() ){
        g1.setExpr(problem["g"].asString());
        g = &g1;
    } else {
        if(problem["Condition Type"].asString() == "mixed"){
            g2.setExpr("x=0", problem["g"]["x=0"][0].asString());
            g2.setExpr("x=1", problem["g"]["x=1"][0].asString());
            g2.setExpr("y=0", problem["g"]["y=0"][0].asString());
            g2.setExpr("y=1", problem["g"]["y=1"][0].asString());
            if(problem["Reigeon Type"].asString()=="Irregular")
                g2.setExpr("D", problem["g"]["D"][0].asString());
            bonDtil.setBondary("x=0", problem["g"]["x=0"][1].asString());
            bonDtil.setBondary("x=1", problem["g"]["x=1"][1].asString());
            bonDtil.setBondary("y=0", problem["g"]["y=0"][1].asString());
            bonDtil.setBondary("y=1", problem["g"]["y=1"][1].asString());
            if(problem["Reigeon Type"].asString()=="Irregular")
                bonDtil.setBondary("D", problem["g"]["D"][1].asString());
        } else {
            g2.setExpr("x=0", problem["g"]["x=0"].asString());
            g2.setExpr("x=1", problem["g"]["x=1"].asString());
            g2.setExpr("y=0", problem["g"]["y=0"].asString());
            g2.setExpr("y=1", problem["g"]["y=1"].asString());
            if(problem["Reigeon Type"].asString()=="Irregular")
                g2.setExpr("D", problem["g"]["D"].asString());
        }
        g2.addConstant("D", "cx", cx);
        g2.addConstant("D", "cy", cy);
        g2.addConstant("D", "R", R);
        g2.setCircle(cx, cy, R);
        g = &g2;
    }

    cout << "----------------------------------------------------------" << endl;
    cout << "Solving..." << endl;

    int cur = clock();
    solver.solve(f, *g, m, problem["Condition Type"].asString(), cx, cy, R, bonDtil);
    double timecost = (double)(clock()-cur)/CLOCKS_PER_SEC;
    timecost = round(timecost*1000)/1000;
    cout << "Solved in " << timecost << "s" << endl;
}

void BVP::output(std::ostream &out){
    cout << "----------------------------------------------------------" << endl;
    out << solver << endl;
    cout << "Results have been saved to result.txt" << endl;
}

void BVP::output(const std::string &ofile){
    ofstream fout(ofile);
    output(fout);
}

void BVP::output(){
    output(cout);
}

std::vector<double> BVP::checkError(){
    if(!problem["Error Check"].asBool())
        return std::vector<double>();
    Fun2D u;
    cout << "----------------------------------------------------------" << endl;
    cout << "u(x,y) = " << problem["u"].asString() << endl;
    u.setExpr(problem["u"].asString());
    if(solver.isPureNeumann()) solver.calcNeumannC(u);

    std::vector<double> err;
    err.push_back(solver.checkError(u, Norm_p(1)));
    cout << "1-norm error: " << err.back() << endl;

    err.push_back(solver.checkError(u, Norm_p(2)));
    cout << "2-norm error: " << err.back() << endl;

    err.push_back(solver.checkError(u, Norm_inf()));
    cout << "max-norm error: " << err.back() << endl;

    return err;
}