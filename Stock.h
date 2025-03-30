// so we need to make a stock price so we can simulate the jumps 


#ifndef Stock.H
#include <vector> 

class Stock{

    protected:
    double sigma;
    double riskfree;

    public:
    Stock(double rf_, double sigma_) : riskfree(rf_), sigma(sigma_) {}

};

class JumpSize{
    
    private:
    double JumpSigma;
    double Jumpmu;
    double k;

    public:
    JumpSize(double mu_ , double sigma_ ) : JumpSigma(simga_), JumpMu(mu_)   { SetK();}
    void SetK();
    double GetK(){return k;}
    

}

class MJD : public Stock{
    private:
    double lambda; 
    double c;
    double k;
    public:
    MJD(double riskfree, double sigma_,double lambda_, double Jumpmu, double JumpSig) : Stock(riskfree,sigma_), lambda(lambda_) , JumpSize(Jumpmu, JumpSig) { SetC(); k = JumpSize.getK();}
    std::vector<double> JumpTimes();     
   
    double ContinuousDynamics(double StartValue, double t1, double t2,);
    double JumpDynamic //this just calls the member function of the class 
    void SetC();

    

};


#endif


