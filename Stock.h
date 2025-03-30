// so we need to make a stock price so we can simulate the jumps 


#ifndef Stock.H
#include <vector> 

class Stock{

    protected:
    double sigma;
    double mu;

    public:
    Stock(double mu_, double sigma_) : mu(mu_), sigma(sigma_) {}

};

class JumpSize{
    
    private:
    double JumpSigma;
    double Jumpmu;
    double k;

    public:
    JumpSize(double mu_ , double sigma_ ) : JumpSigma(simga_), JumpMu(mu_)   { SetK();}
    void SetK();
    

}

class MJD : public Stock{
    private:
    double lambda; 
    double c;
    public:
    MJD(double mu_, double sigma_,double lambda_, double Jumpmu, double JumpSig) : Stock(mu_,sigma_), lambda(lambda_) , JumpSize(Jumpmu, JumpSig) {}
    std::vector<double> JumpTimes();     
   
    double ContinuousDynamics(double StartValue, double t1, double t2,);
    double JumpDynamic //this just calls the member function of the class 
    

    

};


#endif


