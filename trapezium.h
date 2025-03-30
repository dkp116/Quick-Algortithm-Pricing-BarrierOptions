#ifndef TRAPEZIUM_H
#define TRAPEZIUM_H

#define _USE_MATH_DEFINES
#include <cmath>

class Function {
public:
    virtual double evaluate(double t) = 0; // Pure virtual function
    virtual ~Function() = default;
};

class GI : public Function {
private:
    double a, b, H, c, sigma; // Renamed "gamma" to avoid conflict with function
    double gamma(double T1, double T2); // Private helper function

public:
    GI(double a_, double b_, double H_, double c_, double sigma_);
    double evaluate(double t) override; // Override base class method
    // Add T1 and T2 as parameters if needed (adjust base class if necessary)
    double evaluate_gi(double t, double T1, double T2); // Additional method
};

#endif