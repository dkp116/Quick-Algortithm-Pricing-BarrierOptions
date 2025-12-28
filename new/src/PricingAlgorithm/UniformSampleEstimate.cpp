#include "PricingAlgorithm/UniformSampleEstimate.h"
#include <cmath>
#include <random>
#include <iostream>

double UniformSample::NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option, double A,double B, double t1, double t2)
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

double UniformSample::gamma(std::shared_ptr<MertonJumpDynamics> mertonDynamics, double a, double b, double T1, double T2) {
    double c = mertonDynamics->GetC();
    double sigma = mertonDynamics->GetSigma();
    return (1.0 / (std::sqrt(2 * M_PI * (T2 - T1))* sigma)) 
         * std::exp(-(std::pow((a - b) + c * (T2 - T1), 2.0)) 
                     / (2 * sigma * sigma * (T2 - T1)));
}

double UniformSample::evaluate_gi(  std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option,  double a, double b, double t, double T1, double T2) {     //Density of Crossing for the first time during the Brownian Bridge
    double c = mertonDynamics->GetC();    
    double sigma = mertonDynamics->GetSigma();
    double gamma_val = gamma(mertonDynamics,a,b,T1, T2);
    double section1 = ((a - std::log(downAndOut_-> GetBarrier())) / (2 * gamma_val * M_PI * sigma * sigma))
                    * std::pow(t - T1, -3.0/2.0) 
                    * std::pow(T2 - t, -1.0/2.0);

    double expTerm1 = (std::pow((b - std::log(downAndOut_-> GetBarrier()) - c * (T2 - t)), 2.0))
                    / (2 * (T2 - t) * sigma * sigma);

    double expTerm2 = (std::pow((a - std::log(downAndOut_-> GetBarrier()) + c * (t - T1)), 2.0))
                    / (2 * (t - T1) * sigma * sigma);

    return section1 * std::exp(-(expTerm1 + expTerm2));
}

double ExtentionOfTimeInterval(std::vector<double> jumpTimes, double probabilityOfCrossingWithinInterval, double currentJumpInterval){
    return (jumpTimes[currentJumpInterval+1]- jumpTimes[currentJumpInterval]) / (1.0-probabilityOfCrossingWithinInterval);
}

double UniformSample::OneCycle() {
    std::vector<double> jumpTimes;
    jumpTimes = mertonDynamics_->createJumpTimes();      //generates exponenially distributed jump times
    double StockPriceAfterJump = stock_->GetLogStartPrice();
    int currentJumpInterval = 0;
    double StockPriceBeforeJump;
    while(currentJumpInterval+1 < jumpTimes.size()){
      StockPriceBeforeJump = mertonDynamics_->ContinuousDynamics(StockPriceAfterJump,jumpTimes[currentJumpInterval],jumpTimes[currentJumpInterval+1]); //returns stock value at the end of the continous interval 
      double SizeOfJump = mertonDynamics_->Jumpsize();
      long double probabilityOfCrossingWithinInterval = NoCrossingDensity(mertonDynamics_ , option_ , StockPriceAfterJump, StockPriceBeforeJump,jumpTimes[currentJumpInterval],jumpTimes[currentJumpInterval+1]);
      double extentionOfTimeInterval = ExtentionOfTimeInterval(jumpTimes,probabilityOfCrossingWithinInterval,currentJumpInterval);
      std::uniform_real_distribution <> d{jumpTimes[currentJumpInterval], jumpTimes[currentJumpInterval]+extentionOfTimeInterval}; 
       double Sample = d(RandomGenerator::getGenerator());
   
       if(Sample < jumpTimes[currentJumpInterval+1] )   //if there is a crossing during the bridge
       {
        double Payoff = evaluate_gi(mertonDynamics_ ,option_, StockPriceAfterJump,StockPriceBeforeJump, Sample,jumpTimes[currentJumpInterval],jumpTimes[currentJumpInterval+1] ) 
                            * std::exp(-mertonDynamics_->GetRiskFree() * Sample) * option_->GetRebate() * extentionOfTimeInterval;    
        return Payoff;
       }

    if(currentJumpInterval + 2 < jumpTimes.size()){
         StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump ; 
    }
    
      if(StockPriceAfterJump <= std::log(downAndOut_-> GetBarrier()))    //if there is a crossing during the jump
       { 
        double Payoff = std::exp( - mertonDynamics_->GetRiskFree() * jumpTimes[currentJumpInterval+1]) * option_-> GetRebate();         
        return Payoff; 
       }
        
      currentJumpInterval++;

    }

    double TerminalValue = std::exp(StockPriceBeforeJump);
    return option_-> GetRebate() * std::exp(- mertonDynamics_->GetRiskFree() ) * option_->Payoff(TerminalValue) ;   

}

double UniformSample::Price() {
    double price = 0;

    for( int currentJumpInterval = 0 ; currentJumpInterval < iteration_ ; currentJumpInterval++){
        price += OneCycle();
    }

    return price/iteration_;
}



