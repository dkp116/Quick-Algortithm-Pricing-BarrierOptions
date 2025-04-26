#ifndef EstimateGI_h

#define EstimateGI_h


#include <cmath>

struct ModelParams {
    double r = 0.0;
    double T1 = 0.0;
    double T2 = 1.0;
    double X1 = 0.0;
    double X2 = 0.0;
    double LogBarrier = 0.0;
    double sigma = 0.2;

    double time() const { return T2 - T1; }
};
double normal_cdf(double x);
long double A1(const ModelParams& p);
long double A2(const ModelParams& p);
long double C1(const ModelParams& p);
long double C2(const ModelParams& p);
long double B(const ModelParams& p);
long double EstimateGI(const ModelParams& p);

#endif
