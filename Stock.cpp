//stock.cpp
#include "Stock.h"
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
std::vector <double> MJD::JumpTimes(){

    std::vector<double> Times;   //how do i determine how many rv to generate? 
    //want to generate rv between 0 and 1 to represent waitiing times between each step
    //i guess keep generating until between 0 and 1 and subtract from 1 and if goes below or == 0 we stop

    std::random_device rd;
    std::mt19937 gen(rd());

    std::exponential_distribution<> exp_dis(lambda);

    double count = 0.0;

    while(count<=1){
        double time  = exp_dis(gen);
        count = count + time;
        if(count <=1){
            Times.push_back(time);      //this array will be in backwards order surely 
        }
    }
    return Times;
}


double MJD::ContinuousDynamics(double Start , double t1, double t2){
    //so we need to use the guassian disttribution 
    double mean = Start + 
}

void Jumpsize::SetK(){
    k = std::exp(JumpMu-(JumpSigma*JumpSigma)/2 - 1);
}
