#ifndef MertonJumpDynamics_H
#define MertonJumpDynamics_H
#include "IDynamics.h"
#include "Random_Generator.h"


class MertonJumpDynamics : public IDynamics{
    private:
    double lambda_;  //Frequency of jump 
    double c;
    double expectedValueJump; // Jump size constant
    double sigma_;
    double riskfree_;
    double drift;
    public:
    void SetDrift(){ c = riskfree_ - (sigma_ * sigma_ * 0.5) - (lambda_ * expectedValueJump); }
    MertonJumpDynamics(double riskfree, double sigma,
        double lambda, double Jumpmu, double JumpSig) :  
        riskfree_(riskfree), lambda_(lambda), sigma_(sigma) {};
    
    double evolve(double TimeIncrement) override;
    double Jumpsize(double JumpMu , double JumpSigma);
    std::vector<double> createJumpTimes(); 
};

#endif