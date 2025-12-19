
#include "Dynamics/MertonJumpDynamics.h"

double MertonJumpDynamics::evolve(double TimeIncrement){
    std::normal_distribution<> d{0.0, 1.0};
    double generate = d(RandomGenerator::getGenerator());

    return  std::exp(( riskfree_  - 0.5 * sigma_ * sigma_) * TimeIncrement +
                            sigma_ * std::sqrt(TimeIncrement) * generate);   
    }


 double MertonJumpDynamics::Jumpsize(){      
    std::normal_distribution <> d(jumpMu_, jumpSigma_);       
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


double MertonJumpDynamics::ContinuousDynamics(double Start , double t1, double t2){            
    
    double time = t2 - t1;
    double mean = Start + c_ * time;
    double stddev = sigma_ * std::sqrt(time);
    std::normal_distribution<> d{mean, stddev};
    double generate = d(RandomGenerator::getGenerator());
    return generate;      
}