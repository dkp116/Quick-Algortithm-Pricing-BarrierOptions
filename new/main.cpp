#include "Stock.h"
#include "IPricing.h"
#include "IDynamics.h"
#include "StandardMonteCarlo.h"
#include "BlackScholesDynamics.h"
#include "Barrier.h"
#include "Option.h"
#include "MertonJumpDynamics.h"
#include "UniformSampleEstimate.h"
#include "TaylorApproximation.h"
#include <memory>
#include <iostream>

int main() {
    // Create Merton jump diffusion dynamics
    auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0 , 0.1);



    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0);

        

    // Create Brownian Bridge pricing engine


    // UniformSample pricing(s, b, 10000);
    TaylorApproximation pricing2(s, b, 10000);


    // Compute the price
    // double price = pricing.Price();
    double price2 = pricing2.Price();
    // std::cout << "Option price: " << price << std::endl;
    std::cout << "Option price: " << price2 << std::endl;
    return 0;
}
