#ifndef TaylorApproximation_h
#define TaylorApproximation_h
#include "IPricing.h"
#include "Stock_new.cpp"

class TaylorApproximation : public IPricing{
    public:
    TaylorApproximation(double iterations , NewStock stock) : IPricing(iterations, stock) {}
    double Price() override;        
};

#endif