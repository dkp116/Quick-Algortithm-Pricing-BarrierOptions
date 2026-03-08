#include "Stock/Stock.h"
#include "PricingAlgorithm/IPricing.h"
#include "Dynamics/IDynamics.h"
#include "PricingAlgorithm/StandardMonteCarlo.h"
#include "Dynamics/BlackScholesDynamics.h"
#include "Options/Barrier.h"
#include "Options/Option.h"
#include "Dynamics/MertonJumpDynamics.h"
#include "PricingAlgorithm/UniformSampleEstimate.h"
#include "PricingAlgorithm/TaylorApproximation.h"
#include <memory>
#include <iostream>

int main()
{   //so i want to add varience to be calculated here depending on what is inputted 
    // add to constructer 
    // then depending on if its true of false run the same functions
    //seperate by an if statement shouldn't be too bad to add 
    // Create Merton jump diffusion dynamics
    auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0, 0.1);

    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0);

    // Create Brownian Bridge pricing engine

    StandardMonteCarlo pricing (s,b,1000000, StandardErrorCalculation::Included, Time::NotIncluded);
    double  p = pricing.Price();
    std::cout << p << std::endl;
}
