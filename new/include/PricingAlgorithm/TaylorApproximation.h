#ifndef TaylorApproximation_H
#define TaylorApproximation_H

#include "PricingAlgorithm/IPricing.h"
#include "Stock/Stock.h"
#include "Dynamics/IDynamics.h"
#include "Dynamics/MertonJumpDynamics.h"
#include "Options/Barrier.h"

class TaylorApproximation : public IPricing {
private:
    double iteration_;
    std::shared_ptr<IDynamics> stockDynamics_;
    std::shared_ptr<MertonJumpDynamics> mertonDynamics_;
    std::shared_ptr<DownAndOut> downAndOut_;

public:
    TaylorApproximation(std::shared_ptr<Stock> stock,
                   std::shared_ptr<Option> option,
                   double iteration, StandardErrorCalculation isVarienceIncluded, Time isTimeIncluded, VarianceReduction isVarianceReductionIncluded)
        : IPricing(stock, option, iteration, isVarienceIncluded, isTimeIncluded, isVarianceReductionIncluded),
          stockDynamics_(stock_->GetDynamic()) 
    {
       
        mertonDynamics_ = std::dynamic_pointer_cast<MertonJumpDynamics>(stockDynamics_);
        downAndOut_ = std::dynamic_pointer_cast<DownAndOut>(option_);
    }

    bool isThereAJump(double currentJumpInterval, std::vector<double> &jumpTimesFromZeroToOne);
    double NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option, double A,double B, double t1, double t2);
    double OneCycle() override;
    bool isThereCrossingDuringBridge(double stockPriceBeforeJump);
    bool isThereCrossingAfterJump(double stockPriceAfterJump);
    double TerminalValue(double StockPriceBeforeJump, double Pay, double multiplyPi);
    double PriceWithVarianceReduction() override;

};

#endif
