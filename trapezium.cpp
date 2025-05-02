
#include "trapezium.h"

//so we basically want an intergral estiamtor 
//so we want to use trapezium rule for this funciton 

//so we need the original gi(s) evaluation 
//this could be an extension of the option.cpp then 

//so we need to make a function that takes in the same parameters as the estimategi function and the aim is to just see if the difference between the 2

//so takes in the same parameters then okay 
//how to i estimate an intergral with cod e

#include <iostream>
#include <cmath>
#include "Option.h"
using namespace std;

// Define the function to integrate


// Trapezium Rule function
double trapezium_rule(double Time1, double Time2,  double (*f)(MJD, double, double, double, double, double), MJD stock, double a, double b, double T1, double T2,int n) {
    double h = (Time2 - Time1) / n;
    double sum = 0.5 * (f(stock ,a,  b, Time2 , T1, T2) + f(stock ,a,  b, Time1 , T1, T2));  

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(stock ,a,  b, x , T1, T2);
    }

    return sum * h;
}

// int main() {
//     double a = 0.0;       // Lower limit
//     double b = M_PI;      // Upper limit (π)
//     int n = 1000;         // Number of subintervals

//     double result = trapezium_rule(a, b, n);
//     cout << "Approximate integral: " << result << endl;

//     return 0;
// }



