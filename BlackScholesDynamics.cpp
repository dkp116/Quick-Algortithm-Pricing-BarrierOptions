#include "IDynamics.h"

#include "Random_Generator.h"

class BlackScholesDynamics : public IDynamics{
    private:
    double sigma_;
    double riskfree_;
    public:
    BlackScholesDynamics(double riskfree, double sigma) :  riskfree_(riskfree), sigma_(sigma)  {};

    double evolve(double TimeIncrement) override{
    std::normal_distribution<> d{0.0, 1.0};
    double generate = d(RandomGenerator::getGenerator());

    return  std::exp(( riskfree_  - 0.5 * sigma_ * sigma_) * TimeIncrement +
                            sigma_ * std::sqrt(TimeIncrement) * generate);   
    }

};


//we need the stock price to evolve this, so what do we actually get from this class? We will have different implementation of evolve
//depengins on the class we have.