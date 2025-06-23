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
        //there should be a return here surely?
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


//  double Barrier::PriceByBSM(Stock stock){
//     double T = 1.0;
//     double M = 100.0;
//     // double dt = T/M;
//     double dt = 0.005;
//     bool Checker = true;
//     double ValueOfStock = stock.GetS0();
//     double i = 0;
//     // std::cout << "Value of stock issss: " << ValueOfStock << std::endl;
//     while(i < T){
       
//         ValueOfStock = stock.Dynamics(ValueOfStock,dt);
//         // std::cout << "Value of stock is:::: " << ValueOfStock << std::endl;
//         if(ValueOfStock <= H){
//             //logic for payoff for when the barrier is crossed
//             Checker = false;
//             return Rebate  * std::exp(- (i * stock.GetRF())); // feel like we need to multiply by a probability but not sure what this is 


//         }

//         i = i + dt;
//     }

//     if(Checker){
//         //logic for final payoff (normal call)
//         ValueOfStock = stock.Dynamics(ValueOfStock,dt);
//         return  Payoff(ValueOfStock) * std::exp(-stock.GetRF());
//     }
//  }




void Call::Setd1() {
    double S0 = stock.GetS0();
    double K = Strike;
    double r = stock.GetRF();
    double sigma = stock.GetSigma();
    
    d1 = (log(S0 / K) + (r + (sigma * sigma * 0.5)) * T) 
            / (sigma * sqrt(T));
}



void Call::Setd2() {
    double sigma = stock.GetSigma();
    d2 = d1 - sigma * sqrt(T);
}

double N(double x) {
    return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
}

double Call::ClosedPrice() {
    double S0 = stock.GetS0();
    double K = Strike;
    double r = stock.GetRF();
    std::normal_distribution <> d(0.0,1.0);
    
    return N(d1) * S0 - N(d2) * K * exp(-r * T);
}





double DownAndOut::StandardMonteCarlo(MJD stock){

    //okay so we need to simulate the stock price and the jumps and just add in if statements

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

        // Only apply jump if this is not the last step to maturity
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



// double DownAndOut::StandardMonteCarlo(MJD stock) {
//     // Assume T = 1.0.
//     // Immediate knock-out check at t = 0:
//     double StockPrice = stock.GetS0();
//     if (StockPrice <= H) {
//         // Barrier breached immediately; rebate paid at t=0 (no discount).
//         return Rebate;
//     }

//     // Retrieve jump times vector. We assume JumpTimes() returns times in (0,1),
//     // and that if there are no jumps it effectively gives [0,1] or that you handle 0/1 externally.
//     std::vector<double> Times = stock.JumpTimes();
//     // If JumpTimes does not include 0 or 1, and you need to ensure [0,1], you could:
//     //     Times.push_back(0.0);
//     //     Times.push_back(1.0);
//     //     std::sort(Times.begin(), Times.end());
//     // But per your comment, we assume Times already covers [0,1] when there are no jumps.

//     const int TimeStep = 100;       // subdivisions per interval
//     double r = stock.GetRF();

//     // Loop over each interval [Times[i], Times[i+1]]
//     for (int i = 0; i + 1 < static_cast<int>(Times.size()); ++i) {
//         double t_start = Times[i];
//         double t_end   = Times[i + 1];
//         double interval = t_end - t_start;
//         if (interval <= 0.0) {
//             continue;
//         }
//         double dt = interval / TimeStep;

//         // Discrete diffusion steps
//         for (int z = 0; z < TimeStep; ++z) {
//             // Simulate one small step of size dt:
//             StockPrice = stock.Dynamics(StockPrice, dt);
//             // Time at end of this sub-step:
//             double t = t_start + (z + 1) * dt;
//             if (StockPrice < H) {
//                 // Knocked out at this discrete time t:
//                 return Rebate * std::exp(-r * t);
//             }
//         }

//         // After diffusion sub-steps, we are at time t_end (before jump).
//         // If this t_end < 1.0, assume a jump occurs here:
//         if (t_end < 1.0) {
//             double jumpSize = stock.GetJumpDynamics();
//             StockPrice *= std::exp(jumpSize);
//             if (StockPrice < H) {
//                 // Knocked out exactly at jump time t_end:
//                 return Rebate * std::exp(-r * t_end);
//             }
//         }
//         // If t_end == 1.0, it is maturity; we exit loop and handle payoff below.
//     }

//     // If we finish all intervals without knock-out, pay the (discounted) payoff at T=1:
//     double payoff = Payoff(StockPrice); // e.g. max(StockPrice - K, 0)
//     return payoff * std::exp(-r * 1.0);
//}

