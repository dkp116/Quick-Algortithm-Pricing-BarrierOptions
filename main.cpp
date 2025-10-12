//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>
#include <chrono>
#include "Stock_new.cpp"
#include "BlackScholesDynamics.cpp"
#include "MetonJumpDynamics.cpp"
#include <memory>



// Standard normal cumulative distribution function

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


// Analytic formula for a down-and-out European call option (without rebate)

double down_and_out_call(double S, double K, double H, double T, double r, double sigma) {
    if (S <= H) return 0.0; // Option knocked out already

    double lambdas = (r + 0.5 * sigma * sigma) / (sigma * sigma);
    double x1 = std::log(S / H) / (sigma * std::sqrt(T)) + lambdas * sigma * std::sqrt(T);
    double x2 = std::log(H * H / (S * K)) / (sigma * std::sqrt(T)) + lambdas * sigma * std::sqrt(T);

    double vanilla = black_scholes_call(S, K, T, r, sigma);
    double mirror = std::pow(H / S, 2.0 * lambdas) * black_scholes_call(H * H / S, K, T, r, sigma);

    return vanilla - mirror;

    
}



// Down-and-out call with rebate


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



// Merton Jump-Diffusion model call price via Poisson summation 


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




/*
To build for the first time use the cmake file and make sure you are in the 'build' directory to make the program.
Running the program should give the pricing and standard error for each other the methods, note that the timing function will be of little value as all the algorthms are running at the same time. Comment out accordingly to get accurate timings.

For the branches where we use the trapeizum intergral to check refer to the github branch  'Check-Intergral'

https://github.com/dkp116/Quick-Algortithm-Pricing-BarrierOptions.git

*/



int main(){

    std::shared_ptr<IDynamics> dym = std::make_shared<BlackScholesDynamics>(0.1,0.2);

    double start = 100.01;
    NewStock m(start,dym);


}
