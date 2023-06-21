#ifndef _FUNCTION_FACTORY_
#define _FUNCTION_FACTORY_

#include "IVP.h"

class FunctionFactory{
public:
    using CreateFunctionCallBack = TimeFunction* (*)(double);
private:
    using CallbackMap = std::map<std::string, CreateFunctionCallBack>;
public:
    bool registerFunction(const std::string &ID, CreateFunctionCallBack createFn);
    bool unregisterFunction(const std::string &ID);
    TimeFunction* createFunction(const std::string &ID, const double &arg);
private:
    CallbackMap callbacks_;
private:
    FunctionFactory() = default;
    FunctionFactory(const FunctionFactory&) = default;
    FunctionFactory& operator = (const FunctionFactory&) = default;
    ~FunctionFactory() = default;
public:
    static FunctionFactory& Instance();
};

#endif