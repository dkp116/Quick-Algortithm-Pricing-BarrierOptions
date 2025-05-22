//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>
#include <fstream>
#include <iomanip>



int main(){




    MJD stock(100,0.05,0.25,2.0,0.0,0.0);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
   //so we want to put this in a csv
   std::ofstream csvFile("output.csv");
    csvFile << "Time (Years),Stock Price\n";
    csvFile << std::fixed << std::setprecision(6);

    double currentPrice = stock.GetS0();
    double currentTime = 0.0;

    // Write initial value at time 0
    csvFile << currentTime << "," << currentPrice << "\n";

    // Simulate and write each step

    double T = 1.0/1000;
    while(currentTime < 100 * T) {
        currentPrice = stock.Dynamics(currentPrice, T);
        currentTime += T;
        csvFile << currentTime << "," << currentPrice << "\n";
    }

    currentPrice = currentPrice + 5;

    while(currentTime < 600*T) {
        currentPrice = stock.Dynamics(currentPrice, T);
        currentTime += T;
        csvFile << currentTime << "," << currentPrice << "\n";
    }

    currentPrice = currentPrice + 10;

     while(currentTime < 1) {
        currentPrice = stock.Dynamics(currentPrice, T);
        currentTime += T;
        csvFile << currentTime << "," << currentPrice << "\n";
    }



    csvFile.close();
    std::cout << "CSV generated: stock_prices.csv" << std::endl;
    return 0;







    

    

 
   }
