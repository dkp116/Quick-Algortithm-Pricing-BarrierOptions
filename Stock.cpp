//stock.cpp
#include "Stock.h"
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
std::vector<double> MJD::JumpTimes() {
    std::vector<double> Times;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> exp_dis(lambda);

    double count = 0.0;
    Times.push_back(count);
    while (count <= 1) {
        double time = exp_dis(gen);
        if (count + time > 1) break;  // Stop if exceeding 1
        count += time;
        Times.push_back(count);  // Store cumulative jump time
    }
    Times.push_back(1.0);
    return Times;
}

void MJD::SetC(){
    c = riskfree - (sigma * sigma) / 2 - lambda * ExpectedValueJump; //need to calcualte this k value properly 
}



void JumpSize::SetExpectedValueJump() {

    ExpectedValueJump = std::exp(JumpMu + (JumpSigma * JumpSigma) / 2) - 1.0;
}


double MJD::ContinuousDynamics(double Start , double t1, double t2){
    std::random_device rd;
    std::mt19937 gen(rd());
    double time = t2 - t1;
    double mean = Start + c * time;
    double stddev = sigma * std::sqrt(time);
    std::normal_distribution<> d{mean, stddev};
    double generate = d(gen);
    return generate;      
}

double JumpSize::JumpDynamics(){
    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution <> d(JumpMu,JumpSigma);        //the jumps are log normally distribu

    return d(gen);
}   //so we have a size of the jump. now we need to make a function that generates each tree step by step.

std::vector<double> MJD::StockPrices(std::vector<double>times){
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






// #include "Stock.h"
// #include <iostream>
// #include <map>
// #include <random>
// #include <string>
// #include <vector>

// std::vector<double> MJD::JumpTimes() {
//     std::vector<double> Times;
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::exponential_distribution<> exp_dis(lambda);

//     double count = 0.0;
//     Times.push_back(count);
//     while (count <= 1) {
//         double time = exp_dis(gen);
//         if (count + time > 1) break;  // Stop if exceeding 1
//         count += time;
//         Times.push_back(count);  // Store cumulative jump time
//     }
//     Times.push_back(1.0);
//     return Times;
// }

// void MJD::SetC() {
//     c = riskfree + (sigma * sigma) / 2 - lambda * k; 
// }

// void JumpSize::SetK() {
//     k = std::exp(JumpMu + (JumpSigma * JumpSigma) / 2) - 1;
// }

// double MJD::ContinuousDynamics(double Start , double t1, double t2) {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     double time = t2 - t1;
//     double mean = Start + c * time;
//     double stddev = sigma * std::sqrt(time);
//     std::normal_distribution<> d{mean, stddev};
//     return d(gen);      
// }

// double JumpSize::JumpDynamics() {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::normal_distribution<> d(JumpMu, JumpSigma);
//     return d(gen);
// }

// std::vector<double> MJD::StockPrices(std::vector<double> times) {
//     std::vector<double> Prices;
//     double Price = StartPrice;
//     int jumpIndex = 0;  // Index for tracking jump times

//     for (int i = 0; i < times.size() - 1; i++) { 
//         // continuous dynamics 
//         Price = ContinuousDynamics(Price, times[i], times[i+1]);
        
//         // Check if it's time for a jump
//         if (jumpIndex < jump.JumpTimes().size() && times[i] >= jump.JumpTimes()[jumpIndex]) {
//             Price += jump.JumpDynamics();
//             jumpIndex++;  // Move to the next jump time
//         }

//         Prices.push_back(Price);
//     }

//     return Prices;
// }

// void MJD::ScaledJumpTimes(std::vector<double>& JT, double K) {
//     for (double& time : JT) {
//         time = time * K;
//     }
// }
