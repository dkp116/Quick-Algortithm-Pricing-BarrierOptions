#ifndef BROWNIANBRIDGE_H
#define BROWNIANBRIDGE_H
#include "IPricing.h"
#include "Stock_new.cpp"

class BrownianBridge : public IPricing{
   
    public:
    BrownianBridge(double iterations , NewStock stock) : IPricing(iterations,stock) {}
    double Price() override;

};


#endif