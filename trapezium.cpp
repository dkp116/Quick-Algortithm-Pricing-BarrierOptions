// so we need to make an intergral checker when we insert a fucntion with bounds 
// lets make an object called intergral that calculates our intergral using a certain method 


// double f(double x)
// {
//     return // this is our function 
// }
// double TrapeziumMethod(double a , double b , double n)
// {
//     double h = (b-a)/n;
//     double sum = 0.5 * (f(a)+f(b));
//     for(int i = 1; i< n ; i++)
//     {
//         sum += f(a+i*n);
//     }

//     return sum *h // okay this is the trapezium method 
// }

//okay so I want to make a intergration class

// class Intergration{
//     private:
//     double ub, lb;
//     Function f;
//     public:
//     Intergration(Function f_, double ub_, double lb_ ){ub = ub_ ; lb = lb_ ; f = f_;} //constructor 

// };

#include "trapezium.h"

// Constructor (use initializer list)
GI::GI(double a_, double b_, double H_, double c_, double sigma_)
    : a(a_), b(b_), H(H_), c(c_), sigma(sigma_) {}

// Base class method (adjust if T1/T2 are needed)
double GI::evaluate(double t) {
    // Implement or call evaluate_gi with default T1/T2 if needed
    return 0.0;
}

// Custom method for GI-specific logic
double GI::evaluate_gi(double t, double T1, double T2) {
    double gamma_val = gamma(T1, T2);
    
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
double GI::gamma(double T1, double T2) {
    return 1.0 / (2 * M_PI * (T2 - T1) * sigma) 
         * std::exp(-std::pow((a - b) + c * (T2 - T1), 2) 
                     / (2 * sigma * sigma * (T2 - T1)));
}

//so i have a formula and i want to check that this is equal to an intergral 
//so i need a way to approximate that intergral
//so i have written out gi formula already so then i just need a function to estimate gi 



