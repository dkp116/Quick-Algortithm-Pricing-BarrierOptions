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
    std::vector <double> SimulationRange;
    SimulationRange = {1000, 2000 , 5000 , 10000, 20000, 50000, 100000};
    // double simulation = 1000000;
 
//    auto Start = std::chrono::high_resolution_clock::now();
   for(double simulation  : SimulationRange){
       double totalUniform = 0 ;
    double totalTaylor = 0;
    double totalBSM = 0.0;
    double VarienceBSM = 0;
    double BSMVar = 0.0;
    double RealValue = 9.031;
    double UniformVar = 0;
    double TaylorVar = 0;
   
    for( int q =0 ; q<simulation; q++){

    // std::cout << "----------" << std::endl;
    double OneCycle = Derivative.PriceByMJD_Uniform(stock);
    double two  = Derivative.PriceByMJD_Taylor(stock);
    double three = Derivative.PriceByBSM(stock);
    // double four = varReduction + (Derivative.PriceByBSM(stock) - varReduction);
    //  double four = varReduction + (Derivative.PriceByBSM(stock) - varReduction);
    
    UniformVar = UniformVar +  (RealValue - OneCycle) * (RealValue - OneCycle);
    TaylorVar = TaylorVar +  (RealValue - two) * (RealValue - two);
    // VarienceBSM = VarienceBSM +  (RealValue - three) * (RealValue - three);

    totalUniform = totalUniform + OneCycle;
    totalTaylor = totalTaylor + two;
    // totalBSM = totalBSM + three ;
    
    // BSMVar = BSMVar + four;
        // std::cout << three << std::endl;
        
    }    
    // so now we want to price this option without jumps?
    // auto Finish = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(Finish - Start);

    std::cout << "Price using uniform for amount " << simulation << ": " << totalUniform /simulation << std::endl;
     std::cout << "Std uniform: " << std::sqrt((1/(simulation - 1)) * UniformVar) / std::sqrt(simulation) << std::endl;
    

    std::cout << "Price using Taylor for amount " << simulation << ": " <<  totalTaylor/simulation << std::endl;
    std::cout << "Std Taylor : " << std::sqrt((1/(simulation - 1)) * TaylorVar) / std::sqrt(simulation) << std::endl;

    // std::cout << "Price using BSM : " << totalBSM/simulation << std::endl;
    // std::cout << "Std BSM : " << std::sqrt((1/(simulation - 1)) * VarienceBSM) / std::sqrt(simulation) << std::endl;



    // std::cout << "Price using Var reduction BSM : " << BSMVar/10000.0 << std::endl;

    // std::cout <<"Time Taken: " << duration.count() << std::endl;

   }
    

    

 
   }
