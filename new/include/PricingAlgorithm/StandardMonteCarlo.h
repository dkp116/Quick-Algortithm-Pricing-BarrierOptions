#ifndef STANDARDMONTECARLO_H
#define STANDARDMONTECARLO_H
#include "PricingAlgorithm/IPricing.h"
#include "Stock/Stock.h"
#include "Stock/Stock.h"
#include "Dynamics/IDynamics.h"
#include "Dynamics/MertonJumpDynamics.h"
#include "Options/Barrier.h"


class StandardMonteCarlo : public IPricing{
    private:
    double iteration_;
    std::shared_ptr<IDynamics> stockDynamics_;
    std::shared_ptr<MertonJumpDynamics> mertonDynamics_;
    std::shared_ptr<DownAndOut> downAndOut_;
    public:
    StandardMonteCarlo(std::shared_ptr<Stock> stock , 
                        std::shared_ptr<Option> option , 
                        double iteration) :
                         IPricing(stock, option) , 
                         iteration_(iteration), 
                         stockDynamics_(stock->GetDynamic()) {
        mertonDynamics_ = std::dynamic_pointer_cast<MertonJumpDynamics>(stockDynamics_);
        downAndOut_ = std::dynamic_pointer_cast<DownAndOut>(option_);
    }
    double Price() override;  
    double OneCycle() override; 
   

     
};

#endif