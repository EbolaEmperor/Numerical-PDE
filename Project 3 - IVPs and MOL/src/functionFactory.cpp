#include "IVP.h"
#include "Polynomial.h"
#include <iostream>
#include <fstream>
#include "functionFactory.h"
using namespace std;

bool FunctionFactory::registerFunction(const string &ID, CreateFunctionCallBack createFn){
    if(callbacks_.count(ID)) return false;
    callbacks_[ID] = createFn;
    return true;
}

bool FunctionFactory::unregisterFunction(const string &ID){
    if(!callbacks_.count(ID)) return false;
    callbacks_.erase(ID);
    return true;
}

TimeFunction* FunctionFactory::createFunction(const string &ID, const double &arg){
    if(!callbacks_.count(ID)){
        cerr << "FunctionFactory:: No such Function called '" << ID << "'." << endl;
        return nullptr;
    }
    return callbacks_[ID](arg);
}

FunctionFactory& FunctionFactory::Instance(){
    static FunctionFactory factory;
    return factory;
}