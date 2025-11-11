#ifndef STOCK_H
#define STOCK_H
#include "IDynamics.h"
#include <memory>

class Stock {
private:
    double StartPrice;
    std::shared_ptr<IDynamics> dynamics;

public:
    Stock(double S0, std::shared_ptr<IDynamics> dyn)
        : StartPrice(S0), dynamics(std::move(dyn))  {};

    void setDynamics(std::shared_ptr<IDynamics> dyn){
        dynamics = std::move(dyn);
    }
};

#endif