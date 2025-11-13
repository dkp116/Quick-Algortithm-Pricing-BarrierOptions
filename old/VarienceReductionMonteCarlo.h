#ifndef VarienceReductionMonteCarlo_h
#define VarienceReductionMonteCarlo_h
#include "IPricing.h"
#include "Stock_new.cpp"


class VarienceReductionMonteCarlo : public IPricing{
    public:
    VarienceReductionMonteCarlo(double iterations , NewStock stock) : IPricing(iterations, stock) {}
    double Price() override;        
};

#endif