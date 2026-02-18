#include "PricingAlgorithm/StandardMonteCarlo.h"
#include <vector>

//this is for down and out barrier options

double StandardMonteCarlo::Price(){

    std::vector<double> jumpTimesFromZeroToOne;
    jumpTimesFromZeroToOne = mertonDynamics_->createJumpTimes();      //generates exponenially distributed jump times
    double StockPrice= stock_->GetLogStartPrice();
    int i = 0;
    bool Checker = 1;
    double TimeStep = 100;

  for (int i = 0; i + 1 < jumpTimesFromZeroToOne.size(); ++i) {
        double TimeIncrement = jumpTimesFromZeroToOne[i + 1] - jumpTimesFromZeroToOne[i];
        double dt = TimeIncrement / TimeStep;

        for (int z = 0; z < TimeStep; ++z) {
            StockPrice = mertonDynamics_->evolve(dt);
            double t = jumpTimesFromZeroToOne[i] + z * dt;
            if (StockPrice < std::log(downAndOut_->GetBarrier())) {
                return downAndOut_->GetRebate() * std::exp(-mertonDynamics_->GetRiskFree() * t);
            }
        }

        if (i + 1 < jumpTimesFromZeroToOne.size() - 1) {
           
            double SizeOfJump = mertonDynamics_->Jumpsize();
            StockPrice *= std::exp(SizeOfJump);
           
            if (StockPrice < std::log(downAndOut_->GetBarrier())) {
                return downAndOut_->GetRebate() * std::exp(-mertonDynamics_->GetRiskFree() * jumpTimesFromZeroToOne[i + 1]);
            }
        }
    }

    return downAndOut_->Payoff(StockPrice) * std::exp(-mertonDynamics_->GetRiskFree()) * downAndOut_->GetRebate();
}
