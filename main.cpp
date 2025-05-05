//main.cpp
// #include "trapezium.h"
#include <iostream>
#include "Stock.h"
#include "Option.h"
#include "EstimateGI.h"
#include <iostream>

/*
    Tasks:
        1) Finish and clean up proof
        2) Draw diagrams and insert them in
        3) Take screenshots of the method and insert them into the document(optional)
        4) Get started on final pricing formula 
*/

int main(){




    MJD stock(100,0.05,0.25,2.0,0.0,0.1);                //MJD(double initprice, double riskfree, double sigma_,
                                                        //double lambda_, double Jumpmu, double JumpSig)
                                                        // Stock(initprice, riskfree, sigma_), lambda(lambda_), jump(Jumpmu, JumpSig) { SetC(); k = jump.GetK(); }
    DownAndOut Derivative(85,110,1.0);
                                                    //DownAndOut(double H_, double K_, double R_) : Barrier(H_,  K_, R_) {}

    // double totalUniform = 0 ;
    double totalTaylor = 0;
    double Cycles = 10000;
    for( int q =0 ; q<Cycles; q++){

    // std::cout << "----------" << std::endl;
    // double OneCycle = Derivative.PriceByMJD_Uniform(stock);
    double two  = Derivative.PriceByMJD_Taylor(stock);

    // totalUniform = totalUniform + OneCycle;
    totalTaylor = totalTaylor + two;

  
    // std::cout << OneCycle << std::endl;
    }

    // std::cout << "Price using uniform : " << totalUniform /1000000.0 << std::endl;
    std::cout << "Average Error: " << totalTaylor/Cycles<< std::endl;

    


    

 
   
}
