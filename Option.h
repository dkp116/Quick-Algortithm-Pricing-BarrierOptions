#ifndef Option.h
#define Option_h
#include "Stock.h"

class Option{           //assuming T = 1 can adjust for this later on 
    protected:
    double Strike;
    public:
    Option(double k_) : Strike(k_) {}
    
};


class Barrier : public Option {
    protected:
    double H;
    double Rebate;
    public:
    Barrier(double H_, double K_, double R_) : Option(K_) , H(H_), Rebate(R_) {}
    //virtual double CrossingDensity()=0;
    double PriceByMJD_Uniform(MJD stock);
    double PriceByMJD_Taylor(MJD stock);
    virtual double NoCrossingDensity(MJD stock,double A,double B, double t1, double t2) =0;     
   virtual double evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) =0;
   virtual double gamma(MJD stock,double a, double b, double T1, double T2) =0;
   virtual double Payoff(double FinalVal) =0;
   

};


class DownAndOut : public Barrier{
    public:
        DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}
        double Payoff(double FinalVal) override;
        double evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) override;
        double gamma(MJD stock,double a, double b, double T1, double T2) override;


        double NoCrossingDensity(MJD stock,double A,double B, double t1, double t2);  

     

};





#endif