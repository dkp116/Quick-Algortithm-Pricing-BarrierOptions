
#include "MertonJumpDynamics.h"
double MertonJumpDynamics::evolve(double TimeIncrement){
    std::normal_distribution<> d{0.0, 1.0};
    double generate = d(RandomGenerator::getGenerator());

    return  std::exp(( riskfree_  - 0.5 * sigma_ * sigma_) * TimeIncrement +
                            sigma_ * std::sqrt(TimeIncrement) * generate);   
    }


 double MertonJumpDynamics::Jumpsize(double JumpMu , double JumpSigma){      
    std::normal_distribution <> d(JumpMu,JumpSigma);       
     return  d(RandomGenerator::getGenerator());
    }  

std::vector<double> MertonJumpDynamics::createJumpTimes() {      
    std::vector<double> Times;
    std::exponential_distribution<> exp_dis(lambda_);
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