#ifndef STANDARDMONTECARLO_H
#define STANDARDMONTECARLO_H
#include "PricingAlgorithm/IPricing.h"
#include "Stock/Stock.h"

class StandardMonteCarlo : public IPricing{
    private:
    double iteration_;
    public:
    StandardMonteCarlo(std::shared_ptr<Stock> stock , std::shared_ptr<Option> option , double iteration) : IPricing(stock, option) , iteration_(iteration) {}
    double Price() override;        
};

#endif