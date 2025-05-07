//stock.cpp
#include "Stock.h"
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include "Random_Generator.h"


std::vector<double> MJD::JumpTimes() {      //returns the Jump times that are exponentially distributed between [0,1]
    std::vector<double> Times;

    std::exponential_distribution<> exp_dis(lambda);

    double count = 0.0;
    Times.push_back(count);
    while (count <= 1) {
        double time = exp_dis(RandomGenerator::getGenerator());
        if (count + time > 1) break;  // Stop if exceeding 1
        count += time;
        Times.push_back(count);  // Store cumulative jump time
    }
    Times.push_back(1.0);
    return Times;
}



void MJD::SetC(){
    c = riskfree - (sigma * sigma * 0.5) - (lambda * ExpectedValueJump);    
}



void JumpSize::SetExpectedValueJump() {

    ExpectedValueJump = std::exp(JumpMu + (JumpSigma * JumpSigma) * 0.5) - 1.0;     
}


double MJD::ContinuousDynamics(double Start , double t1, double t2){            
    
    assert(t2 >t1 && "times are incorrect");
    double time = t2 - t1;
    double mean = Start + c * time;
    double stddev = sigma * std::sqrt(time);
    std::normal_distribution<> d{mean, stddev};
    double generate = d(RandomGenerator::getGenerator());
    return generate;      
}

double JumpSize::JumpDynamics(){        //returns the Jump Size using normal distribution
    std::normal_distribution <> d(JumpMu,JumpSigma);       
     return  d(RandomGenerator::getGenerator());
}   


double Stock::Dynamics(double Value, double Increment){
     std::normal_distribution<> d{0.0, 1.0};
    double generate = d(RandomGenerator::getGenerator());
  

    
  return Value * std::exp((riskfree - 0.5 * sigma * sigma )*(Increment) + sigma * std::sqrt(Increment) * generate );


}




// std::vector<double> MJD::StockPrices(std::vector<double>times){     //this function is not needed
//     std::vector<double> Prices;
//    double  Price = StartPrice;
//     int i =0;
//     for(int i =0 ; i < times.size() - 1 ; i++){
//         // continuous dynamics 
//         Price = ContinuousDynamics(Price,times[i],times[i+1]);
//         Prices.push_back(Price);
//         Price = Price + jump.JumpDynamics();
//         Prices.push_back(Price);
    
//     }
    
   
//     return Prices;
// }

//  void MJD::ScaledJumpTimes(std::vector<double>& JT, double K){
//     for( double& time : JT){
//         time = time * K;
//     }
//  }





