



#ifndef STOCK_H
#define STOCK_H

#include <vector>
#include <cmath>


class Stock {       //so for the same stock we can price with black scholes.
protected:
    double sigma;
    double riskfree;
    double StartPrice;

public:
    Stock(double initprice, double rf_, double sigma_)
        : StartPrice(initprice), riskfree(rf_), sigma(sigma_) {}
    virtual ~Stock() {}  
   
    double GetS0(){return StartPrice;}
        double GetSigma() { return sigma; }
    double GetRF() const { return riskfree; }
};

// class BlackScholes : public Stock{
//     public:
//     BlackScholes(double S0 , double riskfree_ , double sigma_) : Stock(S0,riskfree_,sigma) {}
//    //so this is the dynamics of the stock where we will incremenetally calculte the price of the stock and see if it violates the barrier.

// };

class JumpSize {
private:
    double JumpSigma;
    double JumpMu;
    double ExpectedValueJump;

public:
    JumpSize(double mu_, double sigma_)
        : JumpSigma(sigma_), JumpMu(mu_) { SetExpectedValueJump(); }
    
    void SetExpectedValueJump();
    double GetExpectedValueJump() { return ExpectedValueJump; }
    double JumpDynamics();
};

class MJD : public Stock {
private:
    double lambda;  //Frequency of jump 
    double c;
    double ExpectedValueJump; // Jump size constant
    JumpSize jump;


public:
    MJD(double initprice, double riskfree, double sigma_,
        double lambda_, double Jumpmu, double JumpSig)
        : Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { ExpectedValueJump = jump.GetExpectedValueJump(); SetC();  }  
     std::vector<double> JumpTimes();
    double Dynamics(double Value, double Increment);
    double GetLogS0() { return std::log(StartPrice); }
    double GetJumpDynamics() { return jump.JumpDynamics(); }
    double ContinuousDynamics(double Start, double t1, double t2);  //  Returns  the endpoint of the continuous process
    void SetC();
    double GetC() { return c; }
    double GetSigma() { return sigma; }
    double GetRF() const { return riskfree; }
    std::vector<double> StockPrices(std::vector<double> times);
    void VisulizeStock();
};

#endif // STOCK_H


    // void ScaledJumpTimes(std::vector<double>& JT, double K);