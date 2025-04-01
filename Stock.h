// so we need to make a stock price so we can simulate the jumps 


#ifndef Stock.H
#include <vector> 

class Stock{

    protected:
    double sigma;
    double riskfree;
    double StartPrice;

    public:
    Stock(double initprice, double rf_, double sigma_) : StartPrice(initprice), riskfree(rf_), sigma(sigma_) {}

};

class JumpSize{
    
    private:
    double JumpSigma;
    double JumpMu;
    double k;

    public:
    JumpSize(double mu_ , double sigma_ ) : JumpSigma(sigma_), JumpMu(mu_)   { SetK();}
    void SetK();
    double GetK(){return k;}
    double JumpDynamics(); 
 
    

};

class MJD : public Stock{
    private:
    double lambda; 
    double c;
    double k;
    JumpSize jump;
    public:
    MJD(double initprice, double riskfree, double sigma_,double lambda_, double Jumpmu, double JumpSig) : Stock(initprice, riskfree,sigma_), lambda(lambda_) , jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK();}
    std::vector<double> JumpTimes();  
    void ScaledJumpTimes(std::vector<double>& JT, double K); 
   
    double ContinuousDynamics(double StartValue, double t1, double t2); 
    void SetC();
    std::vector<double> StockPrices(std::vector<double>times); //this just gets removed 
    
    

};


#endif


