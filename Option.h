#ifndef Option.h
#define Option_h
#include "Stock.h"

class Option{
    protected:
    double k;
    public:
    Option(double k_) : k(k_) {}
    
};


class Barrier : public Option {
    protected:
    double H;
    public:
    Barrier(double H_, double K_) : Option(K_) , H(H_) {}
    virtual double Payoff(){return 0;}
    virtual double CrossingDensity(){return 0;}
    double PriceByMJD(MJD stock);

};


class DownAndOut : public Barrier{
    double Payoff();
    double CrossingDensity();

};





#endif