//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>
#include <chrono>

//with var reduciton what is the formula for this

int main(){




    MJD stock(100,0.05,0.25,2.0,0.0,0.1);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(85,110,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    Call callopt(110,stock);

    double varReduction = callopt.ClosedPrice();

    // std::cout << callopt.ClosedPrice() << std::endl;

    double simulation = 10000;
    double totalUniform = 0;
    double totalTaylor =0;

    for(int i =0 ; i < simulation ; i++){

        double uniform = Derivative.PriceByMJD_Uniform(stock);
        double taylor = Derivative.PriceByMJD_Taylor(stock);

        totalUniform += uniform;
        totalTaylor += taylor;
    }

    std::cout << totalUniform / simulation << std::endl;
    std::cout << totalTaylor/ simulation << std::endl;

    

 
   }
