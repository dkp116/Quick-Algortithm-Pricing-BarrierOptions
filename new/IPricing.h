#ifndef IPricing_h
#define IPricing_h

#include "IDynamics.h"
#include "Stock.h"
#include "Option.h"


class IPricing       
{    
    protected:
    std::shared_ptr<Stock> stock_;
    std::shared_ptr<Option> option_;
    public:
    IPricing(std::shared_ptr<Stock> stock , std::shared_ptr<Option> option ) : stock_(stock), option_(option) {}
    virtual double  Price() = 0;

};                  


#endif