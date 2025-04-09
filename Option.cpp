//Option.cpp

#include "Option.h"
#include <cmath>

// we need to deifine each function lets copy out the formula for gi to help make life alot easier
//this  is our gi function.
// we need to evaluate it a t though 
//a and b would be the prices at the times we have gotten 

//GET H

double DownAndOut::evaluate_gi(MJD stock, double a, double b, double t, double T1, double T2) {
    double c = stock.GetC();
    double sigma = stock.GetSigma();
    double gamma_val = gamma(stock,T1, T2);
    
    double section1 = (a - std::log(H)) / (2 * gamma_val * M_PI * sigma * sigma) 
                    * std::pow(t - T1, -3.0/2.0) 
                    * std::pow(T2 - t, -1.0/2.0);

    // Fixed: Added operators and parentheses
    double expTerm1 = std::pow((b - std::log(H) - c * (T2 - t)), 2.0)
                    / (2 * (T2 - t) * sigma * sigma);

    double expTerm2 = std::pow((b - std::log(H) - c * (t - T1)), 2.0)
                    / (2 * (t - T1) * sigma * sigma);

    return section1 * std::exp(-(expTerm1 + expTerm2));
}

// Gamma function implementation
double DownAndOut::gamma(MJD stock, double a, double b, double T1, double T2) {
    double c = stock.GetC();
    double sigma = stock.GetSigma();
    return 1.0 / (2 * M_PI * (T2 - T1) * sigma) 
         * std::exp(-std::pow((a - b) + c * (T2 - T1), 2) 
                     / (2 * sigma * sigma * (T2 - T1)));
}


double DownAndOut::NoCrossingDensity(MJD stock, double A,double B, double t1, double t2){
  
    double sigma = stock.GetSigma();

    if(B>std::log(H)){
       double  ExpTerm = (2 * (std::log(H) - A) * (std::log(H)-B))
                        / ((t2 - t1) * sigma * sigma);
        return 1 - std::exp(-(ExpTerm));
    }
    else {
        return 0;
    }
}



/*
peudo : Stock jump times -> stock scale jump size 
        
        for(i = 1 <= size )
        {
            But what do we want out of this. We only want to calcualte this as a one of calculation not all at the same time 
            just need to make an if statement and evaluate on the conditions of the crossing.
        }

        We want a way to calculate gi in one instance 
*/