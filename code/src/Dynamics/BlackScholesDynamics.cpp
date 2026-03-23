#include "Dynamics/BlackScholesDynamics.h"

    double BlackScholesDynamics::evolve(double TimeIncrement){
    std::normal_distribution<> d{0.0, 1.0};
    double generate = d(RandomGenerator::getGenerator());

    return  std::exp(( riskfree_  - 0.5 * sigma_ * sigma_) * TimeIncrement +
                            sigma_ * std::sqrt(TimeIncrement) * generate);   
    }




