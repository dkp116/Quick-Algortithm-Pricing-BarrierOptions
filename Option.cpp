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
   
        return Payoff;
       }

    if(i + 2 < Times.size()){
         StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
    
      if(StockPriceAfterJump <= std::log(H))    //if there is a crossing during the jump
       { 
        double Payoff = std::exp( - stock.GetRF() * Times[i+1]) * Rebate;
        Checker = 0;
         
        return Payoff; 
       }
        
      i++;

    }
    if(Checker){    //if there is no crossing for the entire lifespan of the option
        double TerminalValue = std::exp(StockPriceBeforeJump);
       
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
        
     }
     else if(StockPriceAfterJump <= std::log(H)){   //if there is a crossing during the jump
        Checker = 0;
        
        Pay = Pay + Rebate * std::exp(-stock.GetRF() * Times[i+1]) *multiplyPi;
     }
        i++;

    }

    if( Checker){       //if there is no crossing for the entire lifespan of the option
         double TerminalValue = std::exp(StockPriceBeforeJump);
    
        return  Pay + multiplyPi * Payoff(TerminalValue) * std::exp(-stock.GetRF());
    }

 }





double DownAndOut::StandardMonteCarlo(MJD stock){

   

     std::vector<double> Times;
    Times = stock.JumpTimes();      //generates exponenially distributed jump times
    double StockPrice= stock.GetS0();
    int i = 0;
    bool Checker = 1;
    double TimeStep = 100;

  for (int i = 0; i + 1 < Times.size(); ++i) {
        double TimeIncrement = Times[i + 1] - Times[i];
        double dt = TimeIncrement / TimeStep;

        for (int z = 0; z < TimeStep; ++z) {
            StockPrice = stock.Dynamics(StockPrice, dt);
            double t = Times[i] + z * dt;
            if (StockPrice < H) {
                return Rebate * std::exp(-stock.GetRF() * t);
            }
        }

        if (i + 1 < Times.size() - 1) {
           
            double Jumpsize = stock.GetJumpDynamics();
            StockPrice *= std::exp(Jumpsize);
           
            if (StockPrice < H) {
                return Rebate * std::exp(-stock.GetRF() * Times[i + 1]);
            }
        }
    }

    return Payoff(StockPrice) * std::exp(-stock.GetRF()) * Rebate;
}








std::vector<double>Barrier::UniformVarRedCall(MJD stock){
    std::vector<double> Times;
    Times = stock.JumpTimes();      //generates exponenially distributed jump times
    double StockPriceAfterJump = stock.GetLogS0();
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    double BarrierPay = 0.0;
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
        if(BarrierPay == 0.0){
         BarrierPay = evaluate_gi(stock,StockPriceAfterJump,StockPriceBeforeJump, Sample,Times[i],Times[i+1] ) 
                            * std::exp(-stock.GetRF() * Sample) * Rebate * ExtentionOfInterval; 
        }
        Checker = 0;
        
       }

    if(i + 2 < Times.size()){
         StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
    
      if(StockPriceAfterJump <= std::log(H))    //if there is a crossing during the jump
       { 
        if(BarrierPay == 0.0){

        BarrierPay = std::exp( - stock.GetRF() * Times[i+1]) * Rebate;
         
        }
        Checker = 0;
       }
        
      i++;

    }
    if(Checker){    //if there is no crossing for the entire lifespan of the option
        double TerminalValue = std::exp(StockPriceBeforeJump);
        BarrierPay = Rebate * std::exp(- stock.GetRF() ) * Payoff(TerminalValue);
        double Callpay = Payoff(TerminalValue) * std::exp(-stock.GetRF());
       
       return {BarrierPay,Callpay} ; 
    } 
    else{double TerminalValue = std::exp(StockPriceBeforeJump);
        double Callpay = Payoff(TerminalValue) * std::exp(-stock.GetRF());
        return {BarrierPay,Callpay };}

}




