#ifndef MertonJumpDynamics_H
#define MertonJumpDynamics_H
#include "IDynamics.h"
#include "Random_Generator.h"


class MertonJumpDynamics : public IDynamics{
    private:
    double lambda_;  //Frequency of jump 
    double c_;
    double expectedValueJump_; // Jump size constant
    double sigma_;
    double riskfree_;
    double drift;
    double jumpMu_;
    double jumpSigma_;
    public:
    MertonJumpDynamics(double riskfree, double sigma,
        double lambda, double jumpMu, double jumpSigma) :  
        riskfree_(riskfree), lambda_(lambda), sigma_(sigma), jumpSigma_(jumpSigma) , jumpMu_(jumpMu) {};
    
    void SetDrift(){ c_ = riskfree_ - (sigma_ * sigma_ * 0.5) - (lambda_ * expectedValueJump_); }
    double GetC(){return c_;}
    double GetSigma(){return sigma_;}
    double GetRiskFree(){return riskfree_;}
    double evolve(double TimeIncrement) override;
    double Jumpsize();
    std::vector<double> createJumpTimes(); 
    double ContinuousDynamics(double Start , double t1, double t2);
    
};

#endif