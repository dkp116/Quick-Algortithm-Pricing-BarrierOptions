#include "PricingAlgorithm/StandardMonteCarlo.h"
#include <vector>


double StandardMonteCarlo::OneCycle() {
    std::vector<double> jumpTimesFromZeroToOne = mertonDynamics_->createJumpTimes(); 
    double StockPrice = stock_->GetLogStartPrice(); 
    double T = 1.0;                                
    int timeStep = 30;
    bool isThereBarrierCrossing = false;


    for (size_t jumpIncrement = 0; jumpIncrement < jumpTimesFromZeroToOne.size() - 1; ++jumpIncrement) {
        double PayOffForNoCrossing = SimulateStockPath(jumpTimesFromZeroToOne, jumpIncrement, timeStep, StockPrice, isThereBarrierCrossing);
        if (isThereBarrierCrossing)
            return PayOffForNoCrossing;
    }

    return downAndOut_->Payoff(std::exp(StockPrice)) * std::exp(-mertonDynamics_->GetRiskFree() * T);
}

 double StandardMonteCarlo::SimulateStockPath(std::vector<double> &jumpTimesFromZeroToOne, size_t jumpIncrement, int timeStep, double &StockPrice, bool &isThereBarrierCrossing)
{
   
    double timeIncrementBetweenJumps = jumpTimesFromZeroToOne[jumpIncrement + 1] - jumpTimesFromZeroToOne[jumpIncrement];
    double timeStepIncrement = timeIncrementBetweenJumps / timeStep;

    for (int tStep = 0; tStep < timeStep; ++tStep)
    {
        
        double payoffForContinuousCrossing = ContinuousPath(StockPrice, timeStepIncrement, jumpTimesFromZeroToOne, jumpIncrement, tStep, isThereBarrierCrossing);
        if (isThereBarrierCrossing)
            return payoffForContinuousCrossing;
    }

    if (jumpIncrement < jumpTimesFromZeroToOne.size() - 2)
    {
       
        double payOffForJumpCrossing = JumpPath(StockPrice, jumpTimesFromZeroToOne, jumpIncrement, isThereBarrierCrossing);
        if (isThereBarrierCrossing)
            return payOffForJumpCrossing;
    }
    isThereBarrierCrossing = false;
    return {};
}
double StandardMonteCarlo::JumpPath(double &StockPrice, std::vector<double> &jumpTimesFromZeroToOne, size_t jumpIncrement, bool &isThereBarrierCrossing)
{
  
    double jumpSize = mertonDynamics_->Jumpsize();
    StockPrice += jumpSize;
    if (StockPrice < std::log(downAndOut_->GetBarrier()))
    {
        isThereBarrierCrossing = true;
        return downAndOut_->GetRebate() * std::exp(-mertonDynamics_->GetRiskFree() * jumpTimesFromZeroToOne[jumpIncrement + 1]);
    }
    return {};
}

double StandardMonteCarlo::ContinuousPath(double &StockPrice, double timeStepIncrement, std::vector<double> &jumpTimesFromZeroToOne, size_t jumpIncrement, int tStep, bool &isThereBarrierCrossing)
{
   
    StockPrice = std::log(std::exp(StockPrice) * mertonDynamics_->evolve(timeStepIncrement));
    double t = jumpTimesFromZeroToOne[jumpIncrement] + tStep * timeStepIncrement;
    if (StockPrice < std::log(downAndOut_->GetBarrier()))
    {
        isThereBarrierCrossing = true;
        return downAndOut_->GetRebate() * std::exp(-mertonDynamics_->GetRiskFree() * t);
    }
    return {};
}



