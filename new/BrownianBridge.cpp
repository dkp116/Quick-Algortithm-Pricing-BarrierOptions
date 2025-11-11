#include "BrownianBridge.h"

double NoCrossingDensity(Stock stock, double A,double B, double t1, double t2){  //Probability of stock not crossing in the brownian bridge
  
    double sigma = stock_.dynamics.sigma_GetSigma();
    double tau = t2 - t1;
   

  if (B > std::log(H)) {
        double ExpTerm = (2.0 * (std::log(H) - A) * (std::log(H) - B)) / (tau * sigma * sigma);
        return 1.0 - std::exp(-ExpTerm);
    }
    else {
        return 0.0;
    }
}
double BrownianBridge::Price(){
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


// we have a class Brownian Bridge which is only  for down and out options no other option can use this type of pricer, but what if for example there
//exists a brownian bridge method with up and out and down and in each will need their own implementaiton
//create BB(stock) this will tell us the type of pricing algo we need to implement 
// so just have a bunch of override functions ? 
//Price(option == up and out)
//Price(option == down and in ) this will do right 

//there nust be a cleaner implmentation of this that does not involve a switch in terms of have an excliplit override as there will only be a limited amount of options which 
//will use the same pricing method this override can be thought about later, lets just get the implementation of this done for down and out !!