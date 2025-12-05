#ifndef TRAPEZIUM_H
#define TRAPEZIUM_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <functional> 
#include <iostream>

long double trapezium_rule(double Time1, double Time2,  const std::function<double(double)>& f,int n); 

#endif