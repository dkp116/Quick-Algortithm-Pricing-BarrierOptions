//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>
#include <chrono>

//with var reduciton what is the formula for this




double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
}

// Black-Scholes formula for a European call option
double black_scholes_call(double S, double K, double T, double r, double sigma) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / 
                (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

double down_and_out_call(double S, double K, double H, double T, double r, double sigma) {
    if (S <= H) return 0.0; // Option knocked out already

    double lambdas = (r + 0.5 * sigma * sigma) / (sigma * sigma);
    double x1 = std::log(S / H) / (sigma * std::sqrt(T)) + lambdas * sigma * std::sqrt(T);
    double x2 = std::log(H * H / (S * K)) / (sigma * std::sqrt(T)) + lambdas * sigma * std::sqrt(T);

    double vanilla = black_scholes_call(S, K, T, r, sigma);
    double mirror = std::pow(H / S, 2.0 * lambdas) * black_scholes_call(H * H / S, K, T, r, sigma);

    return vanilla - mirror;

    
}




double down_and_out_call_with_rebate(double S, double K, double H, double R, double T, double r, double sigma) {
    if (S <= H) return R * std::exp(-r * T); // Knocked out immediately, rebate paid

    // Compute down-and-out call price without rebate
    double C_DOC = down_and_out_call(S, K, H, T, r, sigma);

    // Rebate term
    double lambda = (r + 0.5 * sigma * sigma) / (sigma * sigma);
    double zeta = std::log(H / S) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);
    double rebate_term = R * std::exp(-r * T) * std::pow(H / S, 2.0 * lambda) * norm_cdf(zeta);

    return C_DOC + rebate_term;
}






double PriceMJD(MJD stock, int N, double Strike) {
    double price = 0.0;

    // Extract parameters
    double S0 = stock.GetS0();
    double muJ = stock.GetJumpMu();       // Mean of log jump size
    double sigJ = stock.GetJumpSig();     // Stddev of log jump size
    double sigma = stock.GetSigma();      // Diffusion volatility
    double r = stock.GetRF();             // Risk-free rate
    double lambda = stock.GetLamda();     // Jump intensity (expected # jumps per year)
    double T = 1.0;                       // Time to maturity in years

    // Compute kappa = E[Y - 1], where Y = e^Z is the jump multiplier
    double kappa = exp(muJ + 0.5 * sigJ * sigJ) - 1.0;

    // Loop over number of jumps
    for (int n = 0; n < N; ++n) {
        // Adjust volatility and drift for n jumps
        double sigma_n = std::sqrt(sigma * sigma + (n * sigJ * sigJ) / T);
        double r_n = r - lambda * kappa + (n * (muJ + 0.5 * sigJ * sigJ)) / T;

        // Poisson probability of n jumps in time T
        double poisson_prob = exp(-lambda * T) * std::pow(lambda * T, n) / std::tgamma(n + 1.0);

        // Black-Scholes price for adjusted parameters
        price += poisson_prob * black_scholes_call(S0, Strike, T,r_n, sigma_n);
    }

    return price;
}








int main(){




    MJD stock(100,0.05,0.25,2.0,0.0,0.1);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(85,110,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    // Call callopt(110,stock);

    // double test = down_and_out_call_with_rebate(100, 110, 90,1.0,1.0,0.05,0.25);     //down_and_out_call_with_rebate(double S, double K, double H, double R, double T, double r, double sigma)

    // double test2 = PriceMJD(stock,10,110);  //bs_jd_call_price(const double S, const double K, const double r, const double sigma, const double T, const int N, const double m, const double lambda, const double nu) {

    // std::cout << test << std::endl;        
    // std::cout << test2 << std::endl;        

  
   double ControlDeri =  PriceMJD(stock,100,110);
     auto start = std::chrono::high_resolution_clock::now();
   std::vector<double> Simulations ={1000,2000,5000,10000,20000,50000,100000};
   for(double simulation : Simulations){

    // double simulation = 100000;
    double totalUniform = 0;
    double totalTaylor =0;
    double totalMonte = 0;
    double MonteVar =0;
    double TaylorVar =0;
    double XX = 0;
    double UU=0;
    double TT =0;
    
    double TotalControl =0;
    double VarControl = 0;
    double XY= 0 ;

    for(int i =0 ; i < simulation ; i++){

        std::vector<double> uniform = Derivative.UniformVarRedCall(stock);

         UU +=   uniform[0] *uniform[0];
        double Control = uniform[1];
        VarControl+= Control * Control;
        TotalControl += Control;
        totalUniform += uniform[0];
        
        XY += uniform[0] * uniform[1];
        // double Taylor = Derivative.PriceByMJD_Taylor(stock);
        // totalTaylor += Taylor;
        // double Monte = Derivative.StandardMonteCarlo(stock);
        // totalMonte += Monte;
        // MonteVar += Monte * Monte;
        // TaylorVar += Taylor * Taylor;
      
    }
    double EUni = totalUniform / simulation;
    double Econtrol  = TotalControl / simulation;
    double EXY  = XY/simulation;
    double VarY  = VarControl/simulation  - Econtrol  * Econtrol;
    double Cov = EXY - EUni * Econtrol;
    double Beta  = Cov / VarY;
    double VarX = (UU/ simulation) - EUni * EUni;

    // double V1 = MonteVar/simulation  - (totalMonte /simulation)  * (totalMonte /simulation);
    // double V2 = TaylorVar/simulation  - (totalTaylor /simulation)  * (totalTaylor/simulation);




    // std::cout << "Monte Carlo " << simulation << ": " << std::sqrt( ((XX / simulation) - Monte * Monte)/simulation )<< std::endl;
    // std::cout << "Uniform " << simulation << ": " << std::sqrt( ((UU/ simulation) - EUni * EUni)/simulation )<< std::endl;
    std::cout << "Varience Control: " <<simulation <<" " <<  std::sqrt(( VarX - 2.0 * Beta * Cov + Beta * Beta * VarY)/simulation) << std::endl;
    // std::cout << "Varience Monte: " <<std::sqrt(V1/simulation) << std::endl;
    // std::cout << "Varience Taylor: " <<std::sqrt(V2/simulation) << std::endl;

    // std::cout << totalMonte / simulation << std::endl;
    // std::cout << totalUniform / simulation << std::endl;
    std::cout << EUni - Beta * (Econtrol - ControlDeri) << std::endl;
  
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
     std::cout << "Function took " << duration.count() << " milliseconds.\n";       //it would just take as long as we need to complete all paths of the stock so the speed we make up actaully goes away .. fuckkk
   


    //there are so many variables it is hard to keep track of what is going on!!!!
  
    
//lets start off by measuring the varience
 
   }

}