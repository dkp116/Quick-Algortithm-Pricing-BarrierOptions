//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include <iostream>

/*
    Tasks:
        1) Finish and clean up proof
        2) Draw diagrams and insert them in
        3) Take screenshots of the method and insert them into the document(optional)
        4) Get started on final pricing formula 
*/

int main(){




    MJD stock(50,0.05,0.3,8.0,0.0,0.05);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        //: Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(45,55,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    double total = 0 ;
    for( int q =0 ; q<1000; q++){

    std::cout << "----------" << std::endl;
    double OneCycle = Derivative.PriceByMJD_Uniform(stock);

    total = total + OneCycle;

  
        // std::cout << OneCycle << std::endl;
        //so the error we are facing is that the payoff is negative  and not sure why this is the case as the formulas i have checked semmed to be correct as well as the workflow need to spend alot more time on this 
    }    

    std::cout << total /1000.0 << std::endl;

    

 
   }
