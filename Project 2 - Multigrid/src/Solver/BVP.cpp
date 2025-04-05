#include "BVP.h"
#include "solver.h"

#include "Core/mathExpr.h"

#include "IntergridOp/Injection.h"
#include "IntergridOp/FullWeighting.h"
#include "IntergridOp/LinearInterpolation.h"
#include "IntergridOp/QuadraticInterpolation.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <json/json.h>
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

void error(const string &er){
    cerr << "[Error] " << er << endl;
    exit(-1);
}

template <int Dim>
void BVP<Dim>::printProblem(){
    checkParse();
    cout << "----------------------------------------------------------" << endl;
    cout << "BVP with " << problem["Condition Type"].asString() << " Condition" << endl;
    if(problem["Dimension"].asInt()==1){
        cout << "-u''(x) = f(x) = " << problem["f"].asString() << endl;
        if(problem["Condition Type"].asString() == "mixed"){
            // x=0
            if(problem["g"]["x=0"][1].asString()=="Dirichlet")
                cout << "u(0) = " << problem["g"]["x=0"][0].asString() << endl;
            else
                cout << "u'(0) = " << problem["g"]["x=0"][0].asString() << endl;
            // x=1
            if(problem["g"]["x=1"][1].asString()=="Dirichlet")
                cout << "u(1) = " << problem["g"]["x=1"][0].asString() << endl;
            else
                cout << "u'(1) = " << problem["g"]["x=1"][0].asString() << endl;
        } else {
            // x=0
            if(problem["Condition Type"].asString()=="Dirichlet")
                cout << "u(0) = " << problem["g"]["x=0"].asString() << endl;
            else
                cout << "u'(0) = " << problem["g"]["x=0"].asString() << endl;
            // x=1
            if(problem["Condition Type"].asString()=="Dirichlet")
                cout << "u(1) = " << problem["g"]["x=1"].asString() << endl;
            else
                cout << "u'(1) = " << problem["g"]["x=1"].asString() << endl;
        }
    } else {
        cout << "f(x,y) = " << problem["f"].asString() << endl;
        if( !problem["g"].isObject() ){
            cout << "g(x,y) = " << problem["g"].asString() << endl;
        } else {
            if(problem["Condition Type"].asString() == "mixed"){
                cout << "g(x,y) = " << problem["g"]["x=0"][0].asString() << "  if x=0  (" << problem["g"]["x=0"][1].asString() << ")" << endl;
                cout << "         " << problem["g"]["x=1"][0].asString() << "  if x=1  (" << problem["g"]["x=1"][1].asString() << ")" << endl;
                cout << "         " << problem["g"]["Down Bondary"][0].asString() << "  if at down-bondary  (" << problem["g"]["Down Bondary"][1].asString() << ")" << endl;
                cout << "         " << problem["g"]["y=1"][0].asString() << "  if y=1  (" << problem["g"]["y=1"][1].asString() << ")" << endl;
            } else {
                cout << "g(x,y) = " << problem["g"]["x=0"].asString() << "  if x=0" << endl;
                cout << "         " << problem["g"]["x=1"].asString() << "  if x=1" << endl;
                cout << "         " << problem["g"]["Down Bondary"].asString() << "  if at down-bondary" << endl;
                cout << "         " << problem["g"]["y=1"].asString() << "  if y=1" << endl;
            }
        }
    }
    cout << "Grid Size: " << problem["Grid Size"].asInt() << endl;
}

template <int Dim>
void BVP<Dim>::read(const string & fname){
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

template <int Dim>
void BVP<Dim>::read(const char * fname){
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

template <int Dim>
void BVP<Dim>::checkParse(){
    if(problem["Condition Type"].asString() != "Dirichlet" && 
       problem["Condition Type"].asString() != "Neumann" && 
       problem["Condition Type"].asString() != "mixed"){
        error("The 'Condition Type' can only be 'Dirichlet', 'Neumann' or 'mixed'.");
    }
    if(problem["Dimension"].asInt() != 1 && problem["Dimension"].asInt() != 2){
        error("The 'Dimension' can only be 1 or 2.");
    }
    if(problem["Dimension"].asInt() != Dim){
        error("The 'Dimension' does not match the template parameter <Dim>");
    }
    int dim = problem["Dimension"].asInt();
    if(dim == 2 && problem["Reigeon Type"].asString() != "Regular" && problem["Reigeon Type"].asString() != "Irregular"){
        error("The 'Reigeon Type' can only be 'Regular' or 'Irregular'.");
    }
    if(!problem["Grid Size"].isInt()){
        error("The 'Grid Size' should be an integer.");
    }
    int n = problem["Grid Size"].asInt();
    do{
        if(n&1) error("The 'Grid Size' must be a positive power of 2.");
        n >>= 1;
    }while(n>1);
    if(!problem["f"].isString()){
        error("Cannot read the function 'f'.");
    }
    if(dim == 1 && !problem["g"].isObject()){
        error("Cannot read the bondary function 'g'.");
    }
    if(dim == 2 && !problem["g"].isString() && !problem["g"].isObject()){
        error("Cannot read the bondary function 'g'.");
    }
    std::vector<std::string> bons;
    bons.push_back("x=0");
    bons.push_back("x=1");
    bons.push_back("y=1");
    bons.push_back("Down Bondary");
    if(problem["g"].isObject()){
        int upn = (dim==1) ? 2 : 4;
        for(int i = 0; i < upn; i++){
            const std::string & p = bons[i];
            if(!problem["g"].isMember(p)){
                error("The bondary function is not complete. Please provide the bondary conditions at '" + p + "'.");
            }
            if(problem["Condition Type"].asString()!="mixed" && !problem["g"][p].isString()){
                error("Cannot correctly read the bondary condition at '" + p + "'. Please provide it as a String.");
            }
            if(problem["Condition Type"].asString()=="mixed" && (!problem["g"][p].isArray() || !problem["g"][p][0].isString()
            || (problem["g"][p][1].asString()!="Dirichlet" && problem["g"][p][1].asString()!="Neumann"))){
                error("Cannot correctly read the bondary condition at '" + p + "'. Please provide it as an String Array of size 2. Tell me it's expression and type, Dirichlet or Neumann.");
            }
        }
    }
    if(problem["Cycle Type"].asString() != "V" && problem["Cycle Type"].asString() != "FMG"){
        error("The 'Cycle Type' can only be 'V' or 'FMG'.");
    }
    if(problem["Restriction"].asString() != "injection" && problem["Restriction"].asString() != "full operator"){
        error("The 'Restriction' can only be 'injection' or 'full operator'.");
    }
    if(problem["Prolongation"].asString() != "linear" && problem["Prolongation"].asString() != "quadradic"){
        error("The 'Prolongation' can only be 'linear' or 'quadradic'.");
    }
    if(problem.isMember("Max Iter") && (!problem["Max Iter"].isInt() || problem["Max Iter"].asInt()<=0)){
        error("The 'Max Iter' should be a positive integer.");
    }
    if(problem.isMember("eps") && (!problem["eps"].isDouble() || problem["eps"].asDouble()<=0)){
        error("The 'eps' should be a positive real number.");
    }
    if(problem.isMember("Error Check") && !problem["Error Check"].isBool()){
        error("The 'Error Check' must be a boolean.");
    }
    if(problem["Check Error"].asBool() && !problem["u"].isString()){
        error("The function 'u' should be correctly provided to check error.");
    }
}

template <>
void BVP<1>::setFunctions(){
    f.setExpr(problem["f"].asString());
    auto derg = std::make_shared<BFuncExpr<1>>();
    if(problem["Condition Type"].asString() == "mixed"){
        derg->setExpr("x=0", problem["g"]["x=0"][0].asString());
        derg->setExpr("x=1", problem["g"]["x=1"][0].asString());
        bonDtil.setBondary("x=0", problem["g"]["x=0"][1].asString());
        bonDtil.setBondary("x=1", problem["g"]["x=1"][1].asString());
    } else {
        derg->setExpr("x=0", problem["g"]["x=0"].asString());
        derg->setExpr("x=1", problem["g"]["x=1"].asString());
    }
    g = std::make_shared<BFuncExpr<1>>(derg);
}

template <>
void BVP<2>::setFunctions(){
    f.setExpr(problem["f"].asString());
    if(problem["9-stencil"].isBool()){
        f.setDeltaExpr(problem["Delta f"].asString());
    }
    if( !problem["g"].isObject() ){
        auto derg = std::make_shared<FuncExpr<2>>();
        derg->setExpr(problem["g"].asString());
        g = std::make_shared<FuncExpr<2>>(derg);
    } else {
        auto derg = std::make_shared<BFuncExpr<2>>();
        if(problem["Condition Type"].asString() == "mixed"){
            derg->setExpr("x=0", problem["g"]["x=0"][0].asString());
            derg->setExpr("x=1", problem["g"]["x=1"][0].asString());
            derg->setExpr("Down Bondary", problem["g"]["Down Bondary"][0].asString());
            derg->setExpr("y=1", problem["g"]["y=1"][0].asString());
            bonDtil.setBondary("x=0", problem["g"]["x=0"][1].asString());
            bonDtil.setBondary("x=1", problem["g"]["x=1"][1].asString());
            bonDtil.setBondary("Down Bondary", problem["g"]["Down Bondary"][1].asString());
            bonDtil.setBondary("y=1", problem["g"]["y=1"][1].asString());
        } else {
            derg->setExpr("x=0", problem["g"]["x=0"].asString());
            derg->setExpr("x=1", problem["g"]["x=1"].asString());
            derg->setExpr("Down Bondary", problem["g"]["Down Bondary"].asString());
            derg->setExpr("y=1", problem["g"]["y=1"].asString());
        }
        g = std::make_shared<BFuncExpr<2>>(derg);
    }
}

template <int Dim>
void BVP<Dim>::solve(){
    checkParse();
    setFunctions();

    cout << "----------------------------------------------------------" << endl;
    cout << "Method: " << problem["Cycle Type"].asString() << "-Cycle" << endl;
    cout << "Restriction: " << problem["Restriction"].asString() << endl;
    cout << "Prolongation: " << problem["Prolongation"].asString() << endl;
    cout << "Prepareing solver..." << endl;

    int cur = clock();
    IntergridOp<Dim> *restri, *prolong;

    if(problem["Restriction"].asString() == "injection")
        restri = new Injection<Dim>;
    else if(problem["Restriction"].asString() == "full operator")
        restri = new FullWeighting<Dim>;
    
    if(problem["Prolongation"].asString() == "linear")
        prolong = new LinearInterpolation<Dim>;
    else if(problem["Prolongation"].asString() == "quadradic")
        prolong = new QuadradicInterpolation<Dim>;

    solver = new Solver<Dim>(*restri, *prolong);
    
    if(problem["9-stencil"].isBool()){
        solver -> useNineStencil();
    }
    if(problem["Reigeon Type"].asString() == "Irregular")
        solver -> setIrregular(true);
    
    int m = problem["Grid Size"].asInt();
    solver -> init(m, f, *g, problem["Condition Type"].asString(), bonDtil);
    
    solver -> setCycle(problem["Cycle Type"].asString());
    if(problem["eps"].isDouble())
        solver -> setEps(problem["eps"].asDouble());
    if(problem["Max Iter"].isInt())
        solver -> setMaxiter(problem["Max Iter"].asInt());
    solver -> solve();
    double timecost = (double)(clock()-cur)/CLOCKS_PER_SEC;
    timecost = round(timecost*1000)/1000;
    cout << "Solved in " << timecost << "s" << endl;
}

template <int Dim>
void BVP<Dim>::output(std::ostream &out){
    cout << "----------------------------------------------------------" << endl;
    solver -> output(out);
    cout << "Results have been saved to result.txt" << endl;
}

template <int Dim>
void BVP<Dim>::output(const std::string &ofile){
    ofstream fout(ofile);
    output(fout);
}

template <int Dim>
void BVP<Dim>::output(){
    output(cout);
}

template <int Dim>
std::vector<double> BVP<Dim>::checkError(){
    if(!problem["Error Check"].asBool())
        return std::vector<double>();
    FuncExpr<Dim> u;
    cout << "----------------------------------------------------------" << endl;

    cout << "u(x,y) = " << problem["u"].asString() << endl;
    u.setExpr(problem["u"].asString());
    
    if(solver->isPureNeumann()) solver->calcNeumannC(u);

    std::vector<double> err;
    err.push_back(solver->checkError(u, Norm_p(1)));
    cout << "1-norm error: " << err.back() << endl;

    err.push_back(solver->checkError(u, Norm_p(2)));
    cout << "2-norm error: " << err.back() << endl;

    err.push_back(solver->checkError(u, Norm_inf()));
    cout << "max-norm error: " << err.back() << endl;

    return err;
}

template class BVP<1>;
template class BVP<2>;