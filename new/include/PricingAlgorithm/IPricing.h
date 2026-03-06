#ifndef IPricing_h
#define IPricing_h

#include "Dynamics/IDynamics.h"
#include "Stock/Stock.h"
#include "Options/Option.h"
#include <iostream>

enum class VarianceCalculation
{
    Included,
    NotIncluded
};
enum class Time
{
    Included,
    NotIncluded
};

class IPricing
{
protected:
    std::shared_ptr<Stock> stock_;
    std::shared_ptr<Option> option_;
    VarianceCalculation varianceCalculation_;
    double iteration_;

public:
    IPricing(std::shared_ptr<Stock> stock, std::shared_ptr<Option> option, double iteration, VarianceCalculation isVarienceIncluded) : stock_(stock), option_(option), iteration_(iteration), varianceCalculation_(isVarienceIncluded) {}
    virtual double OneCycle() = 0;

    void CalculatePriceWithVariance(double &onGoingAverage, double &onGoingSquareAverage)
    {
        double singleCycle = OneCycle();
        onGoingAverage += singleCycle;
        onGoingSquareAverage += singleCycle * singleCycle;
    }
    double Price()
    {
        if (varianceCalculation_ == VarianceCalculation::NotIncluded)
        {

            double price = 0;
            for (int current_cycle = 0; current_cycle < iteration_; current_cycle++)
            {
                price += OneCycle();
            }

            return price / iteration_;
        }

        else
        {
            double onGoingAverage = 0;
            double onGoingSquareAverage = 0;
            for (int current_cycle = 0; current_cycle < iteration_; current_cycle++)
            {
                CalculatePriceWithVariance(onGoingAverage, onGoingSquareAverage);
            }

            double MonteCarloVariance = (onGoingSquareAverage / iteration_) - pow((onGoingAverage / iteration_), 2.0) ;

            std::cout << "Varience is equal to: " << MonteCarloVariance;

            return onGoingAverage / iteration_;
        }
    }
};

#endif