#ifndef STOCK_H
#define STOCK_H
#include "IDynamics.h"
#include <memory>
#include <cmath>

class Stock {
private:
    double StartPrice;
    std::shared_ptr<IDynamics> dynamics;
    double logStartPrice_;

public:
    Stock(double S0, std::shared_ptr<IDynamics> dyn)
        : StartPrice(S0), dynamics(std::move(dyn)), logStartPrice_(std::log(S0)) {};

    void setDynamics(std::shared_ptr<IDynamics> dyn){
        dynamics = std::move(dyn);
    }

     std::shared_ptr<IDynamics> GetDynamic() const  {return dynamics;}

     double GetLogStartPrice(){return logStartPrice_;}
};

#endif