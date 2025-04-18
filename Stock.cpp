//stock.cpp
#include "Stock.h"
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include "random_generator.h"


std::vector<double> MJD::JumpTimes() {
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
    c = riskfree - sigma * sigma * 0.5 - lambda * ExpectedValueJump; //need to calcualte this k value properly 
}



void JumpSize::SetExpectedValueJump() {

    ExpectedValueJump = std::exp(JumpMu + (JumpSigma * JumpSigma) * 0.5) -1.0;
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

double JumpSize::JumpDynamics(){
    

    std::normal_distribution <> d(JumpMu,JumpSigma);        //the jumps are log normally distribu
    
    return  d(RandomGenerator::getGenerator());
}   //so we have a size of the jump. now we need to make a function that generates each tree step by step.





std::vector<double> MJD::StockPrices(std::vector<double>times){     //this function is not needed
    std::vector<double> Prices;
   double  Price = StartPrice;
    int i =0;
    for(int i =0 ; i < times.size() - 1 ; i++){
        // continuous dynamics 
        Price = ContinuousDynamics(Price,times[i],times[i+1]);
        Prices.push_back(Price);
        Price = Price + jump.JumpDynamics();
        Prices.push_back(Price);
    
    }
    
   
    return Prices;
}

 void MJD::ScaledJumpTimes(std::vector<double>& JT, double K){
    for( double& time : JT){
        time = time * K;
    }
 }





