
#ifndef OPTION_H
#define OPTION_H
#include "Stock.h"
#include "EstimateGI.h"

class Option{           
    protected:
    double Strike;
    public:
    Option(double k_) : Strike(k_) {}
    // virtual double Payoff(double FinalVal) =0;
    
};


class Barrier : public Option {
    protected:
    double H;       //barrier 
    double Rebate;
    public:
    Barrier(double H_, double K_, double R_) : Option(K_) , H(H_), Rebate(R_) {}
    double PriceByMJD_Uniform(MJD stock);   //Uniform sampling pricing implementation
    double PriceByMJD_Taylor(MJD stock);    //Approximating the price by taylor expansion on gi intergral
    virtual long double NoCrossingDensity(MJD stock,double A,double B, double t1, double t2) =0;    //Probability of stock not crossing in the brownian bridge
   virtual double evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) =0;    
   virtual double gamma(MJD stock,double a, double b, double T1, double T2) =0; 
   virtual double Payoff(double FinalVal) =0; // what happens to this then
   virtual  long double Call_trapezium(MJD stock, double a, double b, double T1, double T2)=0;
   double PriceByBSM(Stock stock);
};





class Call : public Option{
private:
    double T = 1;
    double d1;
    double d2;
    Stock stock;
    void Setd1();
    void Setd2();
 public:
    Call(double k, Stock s) : Option(k), stock(s) {
        Setd1();
        Setd2();
    }
    double ClosedPrice();
};



class DownAndOut : public Barrier{
    public:
        DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}
        double Payoff(double FinalVal) override;
        double evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) override;
        double gamma(MJD stock,double a, double b, double T1, double T2) override;
       static double CheckGI(MJD stock, double a, double b, double t, double T1, double T2);
        long double Call_trapezium(MJD stock, double a, double b, double T1, double T2) override;
       long double NoCrossingDensity(MJD stock,double A,double B, double t1, double t2) override;  
        //double VarReductionBSM(Call call);

     

};



#endif