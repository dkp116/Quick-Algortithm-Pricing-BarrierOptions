#include "IDynamics.h"

class MertonJump : public IDynamics{
    private:
    double lambda_;  //Frequency of jump 
    double c;
    double ExpectedValueJump; // Jump size constant
    double sigma_;
    double RiskFree_;
    public:
    MertonJump(double riskfree, double sigma,
        double lambda, double Jumpmu, double JumpSig) :  
        RiskFree_(riskfree), lambda_(lambda), sigma_(sigma) {};
};