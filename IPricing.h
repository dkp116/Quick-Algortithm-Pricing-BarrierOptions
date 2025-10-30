#ifndef IPricing_h
#define IPricing_h

#include "IDynamics.h"
#include "Stock_new.cpp"


class IPricing       
{    
private:

public:

};                  

class BrownianBridge : public IPricing{
    private:
    double iterations_;
    NewStock stock_;
    public:
    BrownianBridge(double iterations , NewStock stock) : iterations_(iterations) , stock_(stock) {}

};

class StandardMonteCarlo : public IPricing{
    private:
    double iterations_;
    NewStock stock_;
    public:
    StandardMonteCarlo(double iterations , NewStock stock) : iterations_(iterations) , stock_(stock){}
    

};

class TaylorApproximation : public IPricing{
    private:
    double iterations_;
    NewStock stock_;
    public:
    TaylorApproximation(double iterations , NewStock stock) : iterations_(iterations) , stock_(stock) {}

};

class VarienceReductionMonteCarlo : public IPricing{
    private:
    double iterations_;
    NewStock stock_;
    public:
    VarienceReductionMonteCarlo(double iterations , NewStock stock) : iterations_(iterations) , stock_(stock) {}

};

#endif