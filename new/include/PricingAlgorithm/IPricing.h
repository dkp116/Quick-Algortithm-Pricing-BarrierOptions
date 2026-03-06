#ifndef IPricing_h
#define IPricing_h

#include "Dynamics/IDynamics.h"
#include "Stock/Stock.h"
#include "Options/Option.h"

 enum class VarianceCalculation
    {
        Included,
        NotIncluded
    };
    enum class Time{
        Included,
        NotIncluded
    };

class IPricing       
{    
    protected:
    std::shared_ptr<Stock> stock_;
    std::shared_ptr<Option> option_;
    VarianceCalculation varianceCalculation_;
    
    public:
    IPricing(std::shared_ptr<Stock> stock , std::shared_ptr<Option> option, VarianceCalculation isVarienceIncluded ) : stock_(stock), option_(option), varianceCalculation_(isVarienceIncluded){}
    virtual double  OneCycle() = 0;
    virtual double  Price() = 0;

};                  


#endif