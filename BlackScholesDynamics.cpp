#include "IDynamics.h"

class BlackScholesDynamics : public IDynamics{
    private:
    double sigma_;
    double riskfree_;
    public:
    BlackScholesDynamics(double riskfree, double sigma) :  riskfree_(riskfree), sigma_(sigma)  {};

};