#ifndef TaylorApproximation_H
#define TaylorApproximation_H

#include "IPricing.h"
#include "Stock.h"
#include "IDynamics.h"
#include "MertonJumpDynamics.h"
#include "Barrier.h"

class TaylorApproximation : public IPricing {
private:
    double iteration_;
    std::shared_ptr<IDynamics> stockDynamics_;
    std::shared_ptr<MertonJumpDynamics> mertonDynamics_;
    std::shared_ptr<DownAndOut> downAndOut_;

public:
    TaylorApproximation(std::shared_ptr<Stock> stock,
                   std::shared_ptr<Option> option,
                   double iteration)
        : IPricing(stock, option),
          iteration_(iteration),
          stockDynamics_(stock_->GetDynamic()) 
    {
       
        mertonDynamics_ = std::dynamic_pointer_cast<MertonJumpDynamics>(stockDynamics_);
        downAndOut_ = std::dynamic_pointer_cast<DownAndOut>(option_);
    }


    double NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option, double A,double B, double t1, double t2);
    double OneCycle() override;
    double Price() override;
};

#endif
