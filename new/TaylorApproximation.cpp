 
#include "TaylorApproximation.h"
#include <cmath>
#include <random>
#include <iostream>
#include "EstimateGI.h"
 

double TaylorApproximation::NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option, double A,double B, double t1, double t2)
{  //Probability of stock not crossing in the brownian bridge
  
    double sigma = mertonDynamics->GetSigma();
    double tau = t2 - t1;
   

  if (B > std::log(downAndOut_-> GetBarrier())) {
        double ExpTerm = (2.0 * (std::log(downAndOut_-> GetBarrier()) - A) * (std::log(downAndOut_-> GetBarrier()) - B)) / (tau * sigma * sigma);
        return 1.0 - std::exp(-ExpTerm);
    }
    else {
        return 0.0;
    }
}

 double TaylorApproximation::OneCycle(){
   
    std::vector<double> Times;
    Times = mertonDynamics_->createJumpTimes();   //generates exponenially distributed jump times
    double Pay = 0;
    ModelParams p;
    p.r = mertonDynamics_->GetRiskFree();
    p.sigma = mertonDynamics_ ->GetSigma();
    p.LogBarrier = std::log(downAndOut_->GetBarrier());
    double StockPriceAfterJump = stock_->GetLogStartPrice();
    int i = 0;
    bool Checker = 1;
    double StockPriceBeforeJump = 0.0;
    double multiplyPi = 1;
    while(i+1 < Times.size()){
        StockPriceBeforeJump = mertonDynamics_->ContinuousDynamics(StockPriceAfterJump,Times[i],Times[i+1]);   //returns stock value at the end of the continous interval 
        double SizeOfJump = mertonDynamics_->Jumpsize();    
       long double P_i = NoCrossingDensity(mertonDynamics_ ,option_, StockPriceAfterJump, StockPriceBeforeJump,Times[i],Times[i+1] );    //Probability that there is no corssing during the brownian bridge
        p.T1 = Times[i];
        p.T2 = Times[i+1];
        p.X1 = StockPriceAfterJump;
        p.X2 = StockPriceBeforeJump;
         double J =  EstimateGI(p);
        
        Pay = Pay + option_->GetRebate() * J * multiplyPi;
     if(i + 2 < Times.size()){
        StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
     multiplyPi = multiplyPi * P_i;
     if(StockPriceBeforeJump <= std::log(downAndOut_->GetBarrier())){       //if there is a crossing during the bridge
        Checker = 0;
        return Pay;
        
     }
     else if(StockPriceAfterJump <= std::log(downAndOut_->GetBarrier())){   //if there is a crossing during the jump
        Checker = 0;
        
        return Pay = Pay + option_->GetRebate() * std::exp(- mertonDynamics_->GetRiskFree() * Times[i+1]) *multiplyPi;
     }
        i++;

    }

    if( Checker){       //if there is no crossing for the entire lifespan of the option
         double TerminalValue = std::exp(StockPriceBeforeJump);
    
        return  Pay + multiplyPi * downAndOut_->Payoff(TerminalValue) * std::exp(- mertonDynamics_->GetRiskFree());
    }

 }

  double TaylorApproximation::Price() {

    double price = 0;
    for( int z = 0 ; z < iteration_ ; z++){
        price += OneCycle();
    }

    return price/iteration_;
}

