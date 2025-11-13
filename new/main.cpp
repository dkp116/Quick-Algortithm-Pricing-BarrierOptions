#include "Stock.h"
#include "IPricing.h"
#include "IDynamics.h"
#include "StandardMonteCarlo.h"
#include "BlackScholesDynamics.h"
#include "Barrier.h"
#include "Option.h"
#include "MertonJumpDynamics.h"
#include "UniformSampleEstimate.h"
#include <memory>
#include <iostream>

int main() {
    // Create Merton jump diffusion dynamics
    auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.3, 8, 0 , 0.05);



    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(50.0, dynamic);

    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 55, 45, 1.0);

        

    // Create Brownian Bridge pricing engine
    UniformSample pricing(s, b, 100000);

    // Compute the price
    double price = pricing.Price();
    std::cout << "Option price: " << price << std::endl;
// 
    return 0;
}
