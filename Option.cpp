//Option.cpp

#include "Option.h"
#include <cmath>
#include <random>
#include <iostream>
#include <cassert>
#include <iomanip>
#include "Random_Generator.h"
#include "trapezium.h"


double DownAndOut::evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) {     //Density of Crossing for the first time during the Brownian Bridge
    double c = stock.GetC();    
    double sigma = stock.GetSigma();
    double gamma_val = gamma(stock,a,b,T1, T2);
    double section1 = ((a - std::log(H)) / (2 * gamma_val * M_PI * sigma * sigma))
                    * std::pow(t - T1, -3.0/2.0) 
                    * std::pow(T2 - t, -1.0/2.0);

    double expTerm1 = (std::pow((b - std::log(H) - c * (T2 - t)), 2.0))
                    / (2 * (T2 - t) * sigma * sigma);

    double expTerm2 = (std::pow((a - std::log(H) + c * (t - T1)), 2.0))
                    / (2 * (t - T1) * sigma * sigma);

    return section1 * std::exp(-(expTerm1 + expTerm2));
}

double DownAndOut::gamma(MJD stock, double a, double b, double T1, double T2) {
    double c = stock.GetC();
    double sigma = stock.GetSigma();
    return (1.0 / (std::sqrt(2 * M_PI * (T2 - T1))* sigma)) 
         * std::exp(-(std::pow((a - b) + c * (T2 - T1), 2.0)) 
                     / (2 * sigma * sigma * (T2 - T1)));
}



long double DownAndOut::Call_trapezium(MJD stock, double a, double b, double T1, double T2) {
    auto function = [this, stock, a, b, T1, T2](double t) {
        double r = stock.GetRF();
        return this->evaluate_gi(stock, a, b, t, T1, T2) *  std::exp(- r * t); 
    };
    return trapezium_rule(T1, T2, function, 100);
}



long double DownAndOut::NoCrossingDensity(MJD stock, double A,double B, double t1, double t2){  //Probability of stock not crossing in the brownian bridge
  
    double sigma = stock.GetSigma();
    double tau = t2 - t1;
   

  if (B > std::log(H)) {
        double ExpTerm = (2.0 * (std::log(H) - A) * (std::log(H) - B)) / (tau * sigma * sigma);
        return 1.0 - std::exp(-ExpTerm);
    }
    else {
        return 0.0;
    }
}




double DownAndOut::Payoff(double FinalVal){     
    if(FinalVal > Strike){
        return FinalVal - Strike;

    }

    else {
        return 0.0;
    }
}


double Barrier::PriceByMJD_Uniform(MJD stock){     
    std::vector<double> Times;
    Times = stock.JumpTimes();      //generates exponenially distributed jump times
    double StockPriceAfterJump = stock.GetLogS0();
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    while(i+1 < Times.size()){
      StockPriceBeforeJump = stock.ContinuousDynamics(StockPriceAfterJump,Times[i],Times[i+1]); //returns stock value at the end of the continous interval 
      
      double SizeOfJump = stock.GetJumpDynamics();
      long double P_i = NoCrossingDensity(stock , StockPriceAfterJump, StockPriceBeforeJump,Times[i],Times[i+1] );
      double ExtentionOfInterval = (Times[i+1]- Times[i]) / (1.0-P_i);
      std::uniform_real_distribution <> d{Times[i], Times[i]+ExtentionOfInterval}; 
       double Sample = d(RandomGenerator::getGenerator());
       assert(Sample > Times[i] && "Invalid time of sample ");
   
       if(Sample < Times[i+1] )   //if there is a crossing during the bridge
       {
        double Payoff = evaluate_gi(stock,StockPriceAfterJump,StockPriceBeforeJump, Sample,Times[i],Times[i+1] ) 
                            * std::exp(-stock.GetRF() * Sample) * Rebate * ExtentionOfInterval; 
        Checker = 0;
        std::cout << "MJD:(1) " << Payoff << std::endl;
        return Payoff;
       }

    if(i + 2 < Times.size()){
         StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
    
      if(StockPriceAfterJump <= std::log(H))    //if there is a crossing during the jump
       { 
        double Payoff = std::exp( - stock.GetRF() * Times[i+1]) * Rebate;
        Checker = 0;
         std::cout << "MJD:(2) " << Payoff << std::endl;
        return Payoff; 
       }
        
      i++;

    }
    if(Checker){    //if there is no crossing for the entire lifespan of the option
        double TerminalValue = std::exp(StockPriceBeforeJump);
        std::cout << "MJD: " << Rebate * std::exp(- stock.GetRF() ) * Payoff(TerminalValue) << std::endl;
       return Rebate * std::exp(- stock.GetRF() ) * Payoff(TerminalValue) ; 
    }   

}


    

 double Barrier::PriceByMJD_Taylor(MJD stock){
   
     std::vector<double> Times;
    Times = stock.JumpTimes();  //generates exponenially distributed jump times
    double Pay = 0;
    ModelParams p;
    p.r = stock.GetRF();
    p.sigma = stock.GetSigma();
    p.LogBarrier = std::log(H);
    double StockPriceAfterJump = stock.GetLogS0();
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    double multiplyPi = 1;
    while(i+1 < Times.size()){
        StockPriceBeforeJump = stock.ContinuousDynamics(StockPriceAfterJump,Times[i],Times[i+1]);   //returns stock value at the end of the continous interval 
        double SizeOfJump = stock.GetJumpDynamics();    
       long double P_i = NoCrossingDensity(stock , StockPriceAfterJump, StockPriceBeforeJump,Times[i],Times[i+1] );    //Probability that there is no corssing during the brownian bridge
        p.T1 = Times[i];
        p.T2 = Times[i+1];
        p.X1 = StockPriceAfterJump;
        p.X2 = StockPriceBeforeJump;
        double J =  EstimateGI(p);
      //here we are going to check this:
        // long double Intergral_Check = Call_trapezium(stock, StockPriceAfterJump, StockPriceBeforeJump, Times[i], Times[i+1]);

        // std::cout << "Here is the difference: " << std::abs(J-Intergral_Check) << std::endl;
        Pay = Pay + Rebate * J * multiplyPi;
     if(i + 2 < Times.size()){
        StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
     multiplyPi = multiplyPi * P_i;
     if(StockPriceBeforeJump <= std::log(H)){       //if there is a crossing during the bridge
        Checker = 0;
        return Pay;
        //there should be a return here surely?
     }
     else if(StockPriceAfterJump <= std::log(H)){   //if there is a crossing during the jump
        Checker = 0;
        std::cout << "the Taylor : " << Pay << std::endl;
        Pay = Pay + Rebate * std::exp(-stock.GetRF() * Times[i+1]) *multiplyPi;
     }
        i++;

    }

    if( Checker){       //if there is no crossing for the entire lifespan of the option
         double TerminalValue = std::exp(StockPriceBeforeJump);
           std::cout << "the Taylor : " << Pay + multiplyPi * Payoff(TerminalValue) * std::exp(-stock.GetRF())<< std::endl;
        return  Pay + multiplyPi * Payoff(TerminalValue) * std::exp(-stock.GetRF());
    }

 }


 double Barrier::PriceByBSM(Stock stock){
    //okay so here we shall go through the potential payoffs:
    //how many until we will stop this?
    
    double T = 1;
    double M = 1000;
    double dt = T/M;
    bool Checker = 1;
    double ValueOfStock = stock.GetS0();
    double i = T/M;
    // std::cout << "Value of stock issss: " << ValueOfStock << std::endl;
    while(i <= T){
       
        ValueOfStock = stock.Dynamics(ValueOfStock,dt);
        // std::cout << "Value of stock is:::: " << ValueOfStock << std::endl;
        if(ValueOfStock < H){
            //logic for payoff for when the barrier is crossed
            Checker = 0;
            std::cout << "the BSM: " << Rebate  * std::exp(-i * stock.GetRF()) << std::endl;
            return Rebate  * std::exp(-i * stock.GetRF()); // feel like we need to multiply by a probability but not sure what this is 


        }

        i = i + dt;
    }

    if(Checker){
        //logic for final payoff (normal call)
        std::cout << "The BSM: " <<  std::exp(-stock.GetRF()) * Payoff(ValueOfStock) << std::endl;
        return std::exp(-stock.GetRF()) * Payoff(ValueOfStock);
    }
 }
