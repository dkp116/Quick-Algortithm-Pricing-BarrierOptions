#ifndef IPricing_h
#define IPricing_h

#include "IDynamics.h"
#include "Stock_new.cpp"


class IPricing       
{    
    protected:
    double iterations_;
    NewStock stock_;
public:
    IPricing(double iterations , NewStock stock) : iterations_(iterations) , stock_(stock) {}
    virtual double  Price() = 0;

};                  


#endif