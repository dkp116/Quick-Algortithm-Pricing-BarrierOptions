#include "IDynamics.h"

class MertonJumpDynamics : public IDynamics{
    private:
    double lambda_;  //Frequency of jump 
    double c;
    double ExpectedValueJump; // Jump size constant
    double sigma_;
    double RiskFree_;
    public:
    MertonJumpDynamics(double riskfree, double sigma,
        double lambda, double Jumpmu, double JumpSig) :  
        RiskFree_(riskfree), lambda_(lambda), sigma_(sigma) {};
    
    
};