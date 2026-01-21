
#include "PricingAlgorithm/TaylorApproximation.h"
#include <cmath>
#include <random>
#include <iostream>
#include "PricingAlgorithm/EstimateGI.h"

double TaylorApproximation::NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics, std::shared_ptr<Option> option, double A, double B, double t1, double t2)
{ // Probability of stock not crossing in the brownian bridge

    double sigma = mertonDynamics->GetSigma();
    double tau = t2 - t1;

    if (B > std::log(downAndOut_->GetBarrier()))
    {
        double ExpTerm = (2.0 * (std::log(downAndOut_->GetBarrier()) - A) * (std::log(downAndOut_->GetBarrier()) - B)) / (tau * sigma * sigma);
        return 1.0 - std::exp(-ExpTerm);
    }
    else
    {
        return 0.0;
    }
}

bool TaylorApproximation::isThereAJump(double currentJumpInterval, std::vector<double> &jumpTimesFromZeroToOne)
{
    if (currentJumpInterval == jumpTimesFromZeroToOne.size() - 1)
    {
        return 0;
    }
    return 1;
}

bool TaylorApproximation::isThereCrossingDuringBridge(double stockPriceBeforeJump)
{
    if (stockPriceBeforeJump >= std::log(downAndOut_->GetBarrier()))
    {
        return 0;
    }

    return 1;
}

bool TaylorApproximation::isThereCrossingAfterJump(double stockPriceAfterJump)
{
    if (stockPriceAfterJump >= std::log(downAndOut_->GetBarrier()))
    {
        return 0;
    }

    return 1;
}

double TaylorApproximation::TerminalValue(double StockPriceBeforeJump, double Pay, double multiplyPi)
{
    double TerminalValue = std::exp(StockPriceBeforeJump);

    return Pay + multiplyPi * downAndOut_->Payoff(TerminalValue) * std::exp(-mertonDynamics_->GetRiskFree());
}
double TaylorApproximation::OneCycle()
{

    std::vector<double> jumpTimesFromZeroToOne;
    jumpTimesFromZeroToOne = mertonDynamics_->createJumpTimes(); // generates exponenially distributed jump jumpTimesFromZeroToOne
    double Pay = 0;
    ModelParams p(mertonDynamics_->GetRiskFree(), mertonDynamics_->GetSigma(), std::log(downAndOut_->GetBarrier()));
    double StockPriceAfterJump = stock_->GetLogStartPrice();
    double StockPriceBeforeJump;
    double multiplyPi = 1;
    for (int currentJumpInterval = 0; currentJumpInterval + 1 < jumpTimesFromZeroToOne.size(); currentJumpInterval++)
    {

        StockPriceBeforeJump = mertonDynamics_->ContinuousDynamics(StockPriceAfterJump, jumpTimesFromZeroToOne[currentJumpInterval], jumpTimesFromZeroToOne[currentJumpInterval + 1]);                                                     // returns stock value at the end of the continous interval
        double probabilityOfCrossingWithinInterval = NoCrossingDensity(mertonDynamics_, option_, StockPriceAfterJump, StockPriceBeforeJump, jumpTimesFromZeroToOne[currentJumpInterval], jumpTimesFromZeroToOne[currentJumpInterval + 1]); // Probability that there is no corssing during the brownian bridge
        p.setParameters(jumpTimesFromZeroToOne[currentJumpInterval], jumpTimesFromZeroToOne[currentJumpInterval + 1], StockPriceBeforeJump, StockPriceAfterJump);

        double J = EstimateGI(p);

        Pay = Pay + option_->GetRebate() * J * multiplyPi;
        if (isThereCrossingDuringBridge(StockPriceBeforeJump))
        {

            return Pay;
        }
        if (isThereAJump(currentJumpInterval, jumpTimesFromZeroToOne))
        {
            double SizeOfJump = mertonDynamics_->Jumpsize();
            StockPriceAfterJump = StockPriceBeforeJump + SizeOfJump;
        }
        multiplyPi = multiplyPi * probabilityOfCrossingWithinInterval;

        if (isThereCrossingAfterJump(StockPriceAfterJump))
        {

            return Pay = Pay + option_->GetRebate() * std::exp(-mertonDynamics_->GetRiskFree() * jumpTimesFromZeroToOne[currentJumpInterval + 1]) * multiplyPi;
        }
    }
    return TerminalValue(StockPriceBeforeJump, Pay, multiplyPi);
}

double TaylorApproximation::Price()
{

    double price = 0;
    for (int z = 0; z < iteration_; z++)
    {
        price += OneCycle();
    }

    return price / iteration_;
}
