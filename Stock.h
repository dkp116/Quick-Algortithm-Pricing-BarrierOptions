



#ifndef STOCK_H
#define STOCK_H

#include <vector>
#include <cmath>


class Stock {
protected:
    double sigma;
    double riskfree;
    double StartPrice;

public:
    Stock(double initprice, double rf_, double sigma_)
        : StartPrice(initprice), riskfree(rf_), sigma(sigma_) {}
    virtual ~Stock() {}  // Virtual destructor for base class
};

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
    double lambda;
    double c;
    double ExpectedValueJump; // Jump size constant
    JumpSize jump;


public:
    MJD(double initprice, double riskfree, double sigma_,
        double lambda_, double Jumpmu, double JumpSig)
        : Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { ExpectedValueJump = jump.GetExpectedValueJump(); SetC();  }
    
    std::vector<double> JumpTimes();
    double GetLogS0() { return std::log(StartPrice); }
    void ScaledJumpTimes(std::vector<double>& JT, double K);
    double GetJumpDynamics() { return jump.JumpDynamics(); }
    double ContinuousDynamics(double Start, double t1, double t2);
    void SetC();
    double GetC() { return c; }
    double GetSigma() { return sigma; }
    double GetRF() { return riskfree; }
    std::vector<double> StockPrices(std::vector<double> times);
};

#endif // STOCK_H

