#ifndef _IVP_FAC_H_
#define _IVP_FAC_H_

#include "RungeKutta.h"
#include "LMMs.h"
#include <unordered_map>
using std::unordered_map;

enum SOLVER_TYPE{
    AdamsBashforth,
    AdamsMoulton,
    BDF,
    ClassicalRK,
    ESDIRK,
    GaussLegendre,
    AdaptiveGaussLegendre,
    Fehlberg,
    DormandPrince
};


class Factory{
public:
    static Factory& createFactory()
    {
        static Factory object;
        return object;
    }

public:
    Factory(const Factory&) = delete;
    Factory(Factory&&) = delete;
    Factory& operator=(const Factory&) = delete;
    Factory& operator=(Factory&&) = delete;

private:
    Factory() = default;
    ~Factory() = default;

public:
    TimeIntegrator* createTimeIntegrator(SOLVER_TYPE type){
        switch(type){
            case ClassicalRK:
                return new ClassicalRKSolver();
                break;
            case ESDIRK:
                return new ESDIRKSolver();
                break;
            case Fehlberg:
                return new FehlbergSolver();
                break;
            case DormandPrince:
                return new DormandPrinceSolver();
                break;
            default:
                std::cerr << "[Error] Please provide the order of the method." << std::endl;
                return nullptr;
                break;
        }
    }

    TimeIntegrator* createTimeIntegrator(SOLVER_TYPE type, const int &order){
        switch(type){
            case AdamsBashforth:
                return new AdamsBashforthSolver(order);
                break;
            case AdamsMoulton:
                return new AdamsMoultonSolver(order);
                break;
            case BDF:
                return new BDFSolver(order);
                break;
            case GaussLegendre:
                return new GaussLegendreRKSolver(order);
                break;
            case AdaptiveGaussLegendre:
                return new AdaptiveGaussLegendreRKSolver(order);
                break;
            default:
                return createTimeIntegrator(type);
                break;
        }
    }
};

Factory& TimeIntegratorFactory = Factory::createFactory();

#endif