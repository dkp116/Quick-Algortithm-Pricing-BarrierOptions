//Option.cpp

#include "Option.h"
#include <cmath>
#include <random>
#include <iostream>
#include <cassert>
#include <iomanip>
#include "Random_Generator.h"


//happy with gi but can check again, it appears this is correct
double DownAndOut::evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) {
    double c = stock.GetC();
    double sigma = stock.GetSigma();
    double gamma_val = gamma(stock,a,b,T1, T2);
    
    double section1 = (a - std::log(H)) / (2 * gamma_val * M_PI * sigma * sigma) 
                    * std::pow(t - T1, -3.0/2.0) 
                    * std::pow(T2 - t, -1.0/2.0);

    // Fixed: Added operators and 
    double expTerm1 = std::pow((b - std::log(H) - c * (T2 - t)), 2.0)
                    / (2 * (T2 - t) * sigma * sigma);

    double expTerm2 = std::pow((a - std::log(H) + c * (t - T1)), 2.0)
                    / (2 * (t - T1) * sigma * sigma);

    // std::cout << section1 * std::exp(-(expTerm1 + expTerm2)) << std::endl;
    return section1 * std::exp(-(expTerm1 + expTerm2));
}

// Gamma function implementation, happy with this function
double DownAndOut::gamma(MJD stock, double a, double b, double T1, double T2) {
    double c = stock.GetC();
    double sigma = stock.GetSigma();
    return (1.0 / (std::sqrt(2 * M_PI * (T2 - T1))* sigma)) 
         * std::exp(-(std::pow((a - b) + c * (T2 - T1), 2.0)) 
                     / (2 * sigma * sigma * (T2 - T1)));
}





long double DownAndOut::NoCrossingDensity(MJD stock, double A,double B, double t1, double t2){
  
    double sigma = stock.GetSigma();

   
   

   // std::cout << "B = " << B << std::endl; //
    //std::cout << "Barr = " << Barr << std::endl; //
     double tau = t2 - t1;
    //  std::cout << "tau: " << tau << std::endl;

  if (B > std::log(H)) {
        double ExpTerm = (2 * (std::log(H) - A) * (std::log(H) - B)) / (tau * sigma * sigma);
        // std::cout << "Start is equal to: " << A << std::endl;
        // std::cout << "The barrier is " << std::log(H) << std::endl;
        // std::cout << "End is " << B << std::endl;
       // std::cout << 1.0 - std::exp(-ExpTerm) << std::endl;

        return 1.0 - std::exp(-ExpTerm);
    }
    else {
        // std::cout << 0.000002 << std::endl; 
        return 0.0;
    }
}




double DownAndOut::Payoff(double FinalVal){     //lets assume it is a call option then we can chage the structure to for downandout put 
    if(FinalVal > Strike){
        return FinalVal - Strike;

    }

    else {
        return 0.0;
    }
}


double Barrier::PriceByMJD_Uniform(MJD stock){      //CHANGE FORMULA TO USE LOGARITHM OF INITIAL PRICE 
    std::vector<double> Times;
    Times = stock.JumpTimes();
    //stock.ScaledJumpTimes(Times,k);
  
    

    double StockPriceAfterJump = stock.GetLogS0();
    std::cout << " START ::" << StockPriceAfterJump << std::endl;
    //this should use the logarithm them 
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    while(i+1 < Times.size()){
        std::cout << Times[i] << std::endl;
        std::cout << Times[i+1] << std::endl;
      StockPriceBeforeJump = stock.ContinuousDynamics(StockPriceAfterJump,Times[i],Times[i+1]);
      std::cout << "The  Log Stock Price Before the Jump is : " << StockPriceBeforeJump << std::endl;    //
      double SizeOfJump = stock.GetJumpDynamics();
      long double P_i = NoCrossingDensity(stock , StockPriceAfterJump, StockPriceBeforeJump,Times[i],Times[i+1] );
    //   std::cout<< "Size of Jump is : " << SizeOfJump << std::endl;
       
      
    //   std::cout << "Pi: " << std::setprecision(21) << P_i << std::endl; //
      double ExtentionOfInterval = (Times[i+1]- Times[i]) / (1.0-P_i);
     

      //std::cout << "The Uniform dis times are " << Times[i] << std::endl;
      //std::cout<< "AND : " <<  Times[i]+ExtentionOfInterval << std::endl;
      std::uniform_real_distribution <> d{Times[i], Times[i]+ExtentionOfInterval}; 
       double Sample = d(RandomGenerator::getGenerator());
       assert(Sample > Times[i] && "Invalid time of sample ");
    //    std::cout << "Uniform distibtion: " << Sample << std::endl;  //
       //std::cout << "i counter is : " << i << std::endl;
       if(Sample < Times[i+1] ) //if there is a crossing in the interval then ..  
       {
        // Add logic for evaluating gi * R 
        double Payoff = evaluate_gi(stock,StockPriceAfterJump,StockPriceBeforeJump, Sample,Times[i],Times[i+1] ) // not sure what value this spits out 
                            * std::exp(-stock.GetRF() * Sample) * Rebate * ExtentionOfInterval; 
        Checker = 0;
        
        // std::cout << Payoff << std::endl;
        return Payoff;
        
       }
        StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    std::cout << "The  Log Stock Price After the Jump is : " << StockPriceAfterJump << std::endl; 
      if(StockPriceAfterJump <= std::log(H))
       {
         // Add logic for R * exp-rt
        double Payoff = std::exp( - stock.GetRF() * Times[i+1]) * Rebate;
        Checker = 0;
               //clearly the stock does not make it this far. 
        // std::cout << Payoff << std::endl;
        
        return Payoff;

       
       }
        
      i++;

    }

    if(Checker){
        double TerminalValue = std::exp(StockPriceBeforeJump);
        double Payoffz =  Rebate * std::exp(-stock.GetRF() ) * Payoff(TerminalValue) ;
        // std::cout << Payoffz << std::endl;
       return Rebate * std::exp(- stock.GetRF() ) * Payoff(TerminalValue) ; //how do we get the final price? //the final price would be stockprice before jump how can we change this scope of this variable
    }

    

   

}



 double Barrier::PriceByMJD_Taylor(MJD stock){
    //ok so what is involved in this method:: Generate the stock again
     std::vector<double> Times;
    Times = stock.JumpTimes();
    double Pay = 0;
    //stock.ScaledJumpTimes(Times,k);
    ModelParams p;
    p.r = stock.GetRF();
    p.sigma = stock.GetSigma();
    p.LogBarrier = std::log(H);
    
    

    double StockPriceAfterJump = stock.GetLogS0();
    std::cout << " START ::" << StockPriceAfterJump << std::endl;
    //this should use the logarithm them 
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    double multiplyPi = 1;
    while(i+1 < Times.size()){
        std::cout << Times[i] << std::endl;
        std::cout << Times[i+1] << std::endl;
      StockPriceBeforeJump = stock.ContinuousDynamics(StockPriceAfterJump,Times[i],Times[i+1]);
      std::cout << "The  Log Stock Price Before the Jump is : " << StockPriceBeforeJump << std::endl;    //
      double SizeOfJump = stock.GetJumpDynamics();
      long double P_i = NoCrossingDensity(stock , StockPriceAfterJump, StockPriceBeforeJump,Times[i],Times[i+1] );
     
      
      


      p.T1 = Times[i];
      p.T2 = Times[i+1];
      p.X1 = StockPriceAfterJump;
      p.X2 = StockPriceBeforeJump;

     double J =  EstimateGI(p);

     Pay = Pay + Rebate * J * multiplyPi;
     StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump;

     multiplyPi = multiplyPi * P_i;
     if(StockPriceBeforeJump <= std::log(H)){
        Checker = 0;
        return Pay;
     }

     else if(StockPriceAfterJump <= std::log(H)){
        Checker = 0;
        Pay = Pay + Rebate * std::exp(-stock.GetRF() * Times[i+1]) *multiplyPi;
     }


      
        i++;


    }

    if( Checker){
         double TerminalValue = std::exp(StockPriceBeforeJump);
        return  Pay + multiplyPi * Payoff(TerminalValue) * std::exp(-stock.GetRF());
    }



 }
