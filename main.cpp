// #include "trapezium.h"
#include <iostream>
#include "Stock.h"


int main(){


    // GI func(10, 20, 8, 0.3, 0.05);
    // double result = func.evaluate_gi(0.5,0.3,0.7);

    // std::cout << result << std::endl;


    MJD stock(100,0.9,0.9,4,0.9,0.9);

    std::vector<double>Times =  stock.JumpTimes();
    stock.ScaledJumpTimes(Times,10);
   std::vector<double> Prices = stock.StockPrices(Times);

   for(double price : Prices){
    std::cout<< price << std::endl;
   }

//    std::cout << size << std::endl;

//     std::cout << end << std::endl;
}