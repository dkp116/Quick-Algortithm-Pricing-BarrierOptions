#ifndef UNIFORM_SAMPLE_ESTIMATE_H
#define UNIFORM_SAMPLE_ESTIMATE_H

#include "PricingAlgorithm/IPricing.h"
#include "Stock/Stock.h"
#include "Dynamics/IDynamics.h"
#include "Dynamics/MertonJumpDynamics.h"
#include "Options/Barrier.h"
#include <optional>

class UniformSample : public IPricing
{
private:
    double iteration_;
    std::shared_ptr<IDynamics> stockDynamics_;
    std::shared_ptr<MertonJumpDynamics> mertonDynamics_;
    std::shared_ptr<DownAndOut> downAndOut_;

public:
    UniformSample(std::shared_ptr<Stock> stock,
                  std::shared_ptr<Option> option,
                  double iteration, VarianceCalculation isVarienceIncluded,  Time isTimeIncluded)
        : IPricing(stock, option, iteration , isVarienceIncluded, isTimeIncluded),
          stockDynamics_(stock_->GetDynamic())
    {

        mertonDynamics_ = std::dynamic_pointer_cast<MertonJumpDynamics>(stockDynamics_);
        downAndOut_ = std::dynamic_pointer_cast<DownAndOut>(option_);
    }

    bool isThereAJump(double currentJumpInterval, std::vector<double> &jumpTimesFromZeroToOne);
    double gamma(std::shared_ptr<MertonJumpDynamics> mertonDynamics, double a, double b, double T1, double T2);
    double evaluate_gi(std::shared_ptr<MertonJumpDynamics> mertonDynamics, std::shared_ptr<Option> option, double a, double b, double t, double T1, double T2);
    double NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics, std::shared_ptr<Option> option, double A, double B, double t1, double t2);
    std::optional<double> crossingDuringContinuousIntervalChecker(double StockPriceAfterJump, double StockPriceBeforeJump, std::vector<double> jumpTimesFromZeroToOne, double currentJumpInterval);
    std::optional<double> CrossingDuringJump(double StockPriceAfterJump, std::vector<double> &jumpTimesFromZeroToOne, int currentJumpInterval);
    double OneCycle() override;
};

#endif
