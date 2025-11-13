#ifndef UNIFORM_SAMPLE_ESTIMATE_H
#define UNIFORM_SAMPLE_ESTIMATE_H

#include "IPricing.h"
#include "Stock.h"
#include "IDynamics.h"
#include "MertonJumpDynamics.h"
#include "Barrier.h"

class UniformSample : public IPricing {
private:
    double iteration_;
    std::shared_ptr<IDynamics> stockDynamics_;
    std::shared_ptr<MertonJumpDynamics> mertonDynamics_;
    std::shared_ptr<DownAndOut> downAndOut_;

public:
    UniformSample(std::shared_ptr<Stock> stock,
                   std::shared_ptr<Option> option,
                   double iteration)
        : IPricing(stock, option),
          iteration_(iteration),
          stockDynamics_(stock_->GetDynamic()) 
    {
       
        mertonDynamics_ = std::dynamic_pointer_cast<MertonJumpDynamics>(stockDynamics_);
        downAndOut_ = std::dynamic_pointer_cast<DownAndOut>(option_);
    }

    double gamma(std::shared_ptr<MertonJumpDynamics> mertonDynamics, double a, double b, double T1, double T2);
    double evaluate_gi(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option,  double a, double b, double t, double T1, double T2);
    double NoCrossingDensity(std::shared_ptr<MertonJumpDynamics> mertonDynamics , std::shared_ptr<Option> option, double A,double B, double t1, double t2);
    double OneCycle() override;
    double Price() override;
};

#endif
