#ifndef IPricing_h
#define IPricing_h

#include "Dynamics/IDynamics.h"
#include "Stock/Stock.h"
#include "Options/Option.h"


class IPricing       
{    
    protected:
    std::shared_ptr<Stock> stock_;
    std::shared_ptr<Option> option_;

    enum Varience{
        True,
        False
    };

    enum Time{
        True,
        False
    };
    
    public:
    IPricing(std::shared_ptr<Stock> stock , std::shared_ptr<Option> option ) : stock_(stock), option_(option) {}
    virtual double  OneCycle() = 0;
    virtual double  Price() = 0;

};                  


#endif