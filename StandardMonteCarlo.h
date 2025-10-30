#ifndef STANDARDMONTECARLO_H
#define STANDARDMONTECARLO_H
#include "IPricing.h"
#include "Stock_new.cpp"

class StandardMonteCarlo : public IPricing{
    public:
    StandardMonteCarlo(double iterations , NewStock stock) : IPricing(iterations, stock) {}
    double Price() override;        
};

#endif