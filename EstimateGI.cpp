// //so this is just going to be a formula 
// #include <cmath>


// double normal_cdf(double x) {
//     return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
// }

// long double A1(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     (2 * r)/(sigma) *  (T2 - T1) * (X1 - X2) * std::exp(-(2 * (X1 - LogBarrier)(X2 - LogBarrier)) / ((T2 - T1) * sigma * sigma ));
// }

// long double A2(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     ((2 * r * (T2 - T1)) / (sigma)) * ((X1 - X2) - 2 * LogBarrier);
// }


// long double C1(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     long double expo = (X1 - X2 );
//     std::sqrt(2 * M_PI * (T2 - T1)) * std::exp(pow(expo,2.0) / ( 2 * ( T2 - T1) * sigma * sigma )) 
//                                          * (normal_cdf((X1 + X2 - 2 * LogBarrier)/(std::sqrt(2 * (T2 - T1) * sigma * sigma))) - 1.0 )  ;
// }


// long double C2(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     long double expo = (X1 - X2 );
//     std::sqrt(2 * M_PI * (T2 - T1)) * std::exp(pow(expo,2.0) / ( 2 * ( T2 - T1) * sigma * sigma )) 
//                                         * (normal_cdf((X1 - X2 )/(std::sqrt(2 * (T2 - T1) * sigma * sigma))) - 1.0);
// }


// long double B(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     8 - (2 * r * (T2 - T1)) + ((2 * r) /(sigma * sigma)) * (X1 - X2) * (X1 + X2 - LogBarrier);
// }



// long double EstimateGI(double r, double T1, double T2,double X1, double X2, double LogBarrier, double sigma){
//     double time = T2-T1;
//     if(X2 > LogBarrier){
//         std::exp(-r * T1) * std::exp(-(2(X1-LogBarrier) * (X2 - LogBarrier))/(time * sigma * sigma)) 
//                                     + ((r * (X1-LogBarrier) * std::exp(r * T2 - 2 * r *T1 ))/(8 * sigma))
//                                     * ( A1(r, T1, T2, X1,  X2,  LogBarrier, sigma) + C1() * B());
//     }
//     else{
//         std::exp(-r * T1) + ((r * ( X1 - LogBarrier ) 
//                     * std::exp(r * T2 - 2 * r * T1))/ (8 * sigma) ) 
//                         * ( A2() + C2() * B());
//     }


// }



#include "EstimateGI.h"         //Taylor expnansion estimation of the intergral disounted density of gi 
double normal_cdf(double x) {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

long double A1(const ModelParams& p) {
    return (2 * p.r / p.sigma) * p.time() * (p.X1 - p.X2) *
           std::exp(-2 * (p.X1 - p.LogBarrier) * (p.X2 - p.LogBarrier) /
                    (p.time() * p.sigma * p.sigma));
}

long double A2(const ModelParams& p) {
    return (2 * p.r * p.time()) / p.sigma * ((p.X1 - p.X2) - 2 * p.LogBarrier);
}

long double C1(const ModelParams& p) {
    double expo = (p.X1 - p.X2);
    return std::sqrt(2 * M_PI * p.time()) *
           std::exp(expo * expo / (2 * p.time() * p.sigma * p.sigma)) *
           (normal_cdf((p.X1 + p.X2 - 2 * p.LogBarrier) /
                       std::sqrt(2 * p.time() * p.sigma * p.sigma)) -
            1.0);
}

long double C2(const ModelParams& p) {
    double expo = (p.X1 - p.X2);
    return std::sqrt(2 * M_PI * p.time()) *
           std::exp(expo * expo / (2 * p.time() * p.sigma * p.sigma)) *
           (normal_cdf((p.X1 - p.X2) /
                       std::sqrt(2 * p.time() * p.sigma * p.sigma)) -
            1.0);
}

long double B(const ModelParams& p) {
    return 8 - (2 * p.r * p.time()) +
           (2 * p.r / (p.sigma * p.sigma)) * (p.X1 - p.X2) *
               (p.X1 + p.X2 - p.LogBarrier);
}

long double EstimateGI(const ModelParams& p) {
    if (p.X2 > p.LogBarrier) {
        return std::exp(-p.r * p.T1) *
                   std::exp(-2 * (p.X1 - p.LogBarrier) *
                            (p.X2 - p.LogBarrier) /
                            (p.time() * p.sigma * p.sigma)) +
               ((p.r * (p.X1 - p.LogBarrier) *
                 std::exp(p.r * p.T2 - 2 * p.r * p.T1)) /
                (8 * p.sigma)) *
                   (A1(p) + C1(p) * B(p));
    } else {
        return std::exp(-p.r * p.T1) +
               ((p.r * (p.X1 - p.LogBarrier) *
                 std::exp(p.r * p.T2 - 2 * p.r * p.T1)) /
                (8 * p.sigma)) *
                   (A2(p) + C2(p) * B(p));
    }
}
