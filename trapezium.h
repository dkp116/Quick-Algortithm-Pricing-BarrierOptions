#ifndef TRAPEZIUM_H
#define TRAPEZIUM_H

#define _USE_MATH_DEFINES
#include <cmath>


double trapezium_rule(double Time1, double Time2,  double (*f)(MJD, double, double, double, double, double), MJD stock, double a, double b, double T1, double T2,int n);

#endif