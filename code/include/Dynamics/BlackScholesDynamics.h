#ifndef BlackScholesDynamics_H
#define BlackScholesDynamics_H

#include "Dynamics/IDynamics.h"
#include "RandomGenerator/Random_Generator.h"

class BlackScholesDynamics : public IDynamics{
    private:
    double sigma_;
    double riskfree_;
    public:
    BlackScholesDynamics(double riskfree, double sigma) :  riskfree_(riskfree), sigma_(sigma)  {};

    double evolve(double TimeIncrement) override;


};


#endif