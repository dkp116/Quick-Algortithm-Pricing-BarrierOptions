//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>



int main(){




    MJD stock(100,0.05,0.25,2.0,0.0,0.1);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(85,110,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    double totalUniform = 0 ;
    double totalTaylor = 0;
    double totalBSM = 0;
    for( int q =0 ; q<100000; q++){

    std::cout << "----------" << std::endl;
    double OneCycle = Derivative.PriceByMJD_Uniform(stock);
    double two  = Derivative.PriceByMJD_Taylor(stock);
    // double three = Derivative.PriceByBSM(stock);

    totalUniform = totalUniform + OneCycle;
    totalTaylor = totalTaylor + two;
    // totalBSM = totalBSM + three;
  
        // std::cout << OneCycle << std::endl;
        
    }    
    //so now we want to price this option without jumps?

    std::cout << "Price using uniform : " << totalUniform /100000.0 << std::endl;
    std::cout << "Price using Taylor : " << totalTaylor/100000.0 << std::endl;
    // std::cout << "Price using BSM : " << totalBSM/10000.0 << std::endl;


    

    

 
   }
