// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include <iostream>


int main(){


    // GI func(10, 20, 8, 0.3, 0.05);
    // double result = func.evaluate_gi(0.5,0.3,0.7);

    // std::cout << result << std::endl;


    MJD stock(0.1,0.05,200);

    std::vector<double> Timings = stock.JumpTimes();

    for(double time : Timings){
        std::cout << time << std::endl;
    }
}