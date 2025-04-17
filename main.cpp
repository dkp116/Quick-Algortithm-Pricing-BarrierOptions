// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include <iostream>


int main(){




    MJD stock(100,0.05,0.005,100.0,0.8,0.5);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(95,110,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    double total = 0 ;
    for( int q =0 ; q<10000 ; q++){

    
    double OneCycle = Derivative.PriceByMJD_Uniform(stock);

    total = total + OneCycle;

    // std::cout <<  OneCycle << std::endl;
   

    }    

    std::cout <<  total/10000.0  << std::endl;

 
   }
