
#include "trapezium.h"

using namespace std;

  long double trapezium_rule(double Time1, double Time2, const std::function<double(double)>& f, int n) { //we are inserting the gi funciton
    double h = (Time2 - Time1) / n;
    long double sum = 0.0;

    for (int i = 1; i < n; ++i) {  // We are excluding the endpoint in the intergral approximation as the gi function can not be evaluated at these points 
        double t = Time1 + i * h;
        sum += f(t);
    }

    return sum * h;
}




