#include "IDynamics.h"
#include <memory>

class NewStock {
private:
    double StartPrice;
    std::shared_ptr<IDynamics> dynamics;

public:
    NewStock(double S0, std::shared_ptr<IDynamics> dyn)
        : StartPrice(S0), dynamics(std::move(dyn))  {};

    void setDynamics(std::shared_ptr<IDynamics> dyn){
        dynamics = std::move(dyn);
    }
};

/*
refactoring to do :
1) make the different dynamic classes:
2) allow the stock to take in the diff dynamics, would this impliclty convert i believe so.
3) add the correct functionalility to the dyanmic or stock
*/