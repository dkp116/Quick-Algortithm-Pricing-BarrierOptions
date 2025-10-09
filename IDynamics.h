#ifndef IDYNAMICS_H
#define IDYNAMICS_H
#include "Stock.h"

class IDynamics{
    public:
    virtual ~IDynamics() = default;
    virtual double evolve(Stock stock) =0 ;
};


#endif