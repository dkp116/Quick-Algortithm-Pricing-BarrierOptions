#include "TestDownAndOutPricing.h"
#include "Options/Barrier.h"

double price_down_and_out_call_with_taylor_series() {
    // Use Merton jump diffusion dynamics
    auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0, 0.1);

    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    // Define a Down-and-Out European call option
    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0
    );

    // Create Monte Carlo pricing engine
    TaylorApproximation pricing(s, b, 100000, StandardErrorCalculation::NotIncluded, Time::NotIncluded);  

    // Compute the price
    return pricing.Price();
}

double price_down_and_out_call_with_uniform_distribution() {
    // Use Merton jump diffusion dynamics
    auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0, 0.1);

    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    // Define a Down-and-Out European call option
    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0
    );

    // Create Monte Carlo pricing engine
    UniformSample pricing(s, b, 100000, StandardErrorCalculation::NotIncluded, Time::NotIncluded);  

    // Compute the price
    return pricing.Price();
}

double price_down_and_out_call_with_standard_monte_carlo(){
     auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0, 0.1);

    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    // Define a Down-and-Out European call option
    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0
    );

    // Create Monte Carlo pricing engine
    StandardMonteCarlo pricing(s, b, 10000, StandardErrorCalculation::NotIncluded, Time::NotIncluded);  

    // Compute the price
    return pricing.Price();

}



double price_down_and_out_call_with_standard_monte_carlo_and_varience(){
     auto dynamic = std::make_shared<MertonJumpDynamics>(0.05, 0.25, 2, 0, 0.1);

    // Create stock with starting price and dynamics
    auto s = std::make_shared<Stock>(100.0, dynamic);

    // Define a Down-and-Out European call option
    std::shared_ptr<Option> b = std::make_shared<DownAndOut>(
        ExerciseType::European, OptionType::Call, 110, 85, 1.0
    );

    // Create Monte Carlo pricing engine
    StandardMonteCarlo pricing(s, b, 10000, StandardErrorCalculation::Included, Time::NotIncluded);  

    // Compute the price
    return pricing.Price();

}
