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
    virtual double NoCrossingDensity(MJD stock, double t1, double t2);     
   virtual double evaluate_gi(MJD stock, double t, double T1, double T2);
   virtual double gamma(MJD stock,double T1, double T2);
   

};


class DownAndOut : public Barrier{
    double Payoff();
    double evaluate_gi(MJD stock, double t, double T1, double T2);
    double gamma(MJD stock,double T1, double T2);
    double NoCrossingDensity(MJD stock, double t1, double t2);
    
};





#endif