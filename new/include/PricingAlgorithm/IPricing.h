#ifndef IPricing_h
#define IPricing_h

#include "Dynamics/IDynamics.h"
#include "Stock/Stock.h"
#include "Options/Option.h"
#include <iostream>
#include <chrono>

// so I want to  be able to profile this
// I want to generate a time for the calculations
// the pricing algo should look like what then

// the same just with a time right?
enum class StandardErrorCalculation
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
    StandardErrorCalculation includeStandardError_;
    Time includeTime_;
    double iteration_;
    double standard_error_;
    double time_;

public:
    IPricing(std::shared_ptr<Stock> stock, std::shared_ptr<Option> option,
             double iteration, StandardErrorCalculation isStandardErrorIncluded,
             Time isTimeIncluded) : stock_(stock), option_(option), iteration_(iteration), includeStandardError_(isStandardErrorIncluded), includeTime_(isTimeIncluded) {}
    virtual double OneCycle() = 0;

    void CalculatePriceWithVariance(double &onGoingAverage, double &onGoingSquareAverage)
    {
        double singleCycle = OneCycle();
        onGoingAverage += singleCycle;
        onGoingSquareAverage += (singleCycle * singleCycle);
    }
    double Price()
    {
        if (includeTime_ == Time::Included)
        {

            if (includeStandardError_ == StandardErrorCalculation::NotIncluded)
            {
                auto start = std::chrono::high_resolution_clock::now();

                double priceWithoutVarience = PriceWithoutVarience();
                auto end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> elapsed = end - start;
                time_ = elapsed.count();
                return PriceWithoutVarience();
            }

            else
            {
                auto start = std::chrono::high_resolution_clock::now();

                double priceWithoutVarience = PriceWithVarience();
                auto end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> elapsed = end - start;
                time_ = elapsed.count();

                return PriceWithVarience();
            }
        }
        else
        {

            if (includeStandardError_ == StandardErrorCalculation::NotIncluded)
            {

                return PriceWithoutVarience();
            }

            else
            {
                return PriceWithVarience();
            }
        }
    }
    double PriceWithVarience()
    {
        double onGoingAverage = 0;
        double onGoingSquareAverage = 0;
        for (int current_cycle = 0; current_cycle < iteration_; current_cycle++)
        {
            CalculatePriceWithVariance(onGoingAverage, onGoingSquareAverage);
        }

        double variance =
            (onGoingSquareAverage / iteration_) -
            pow(onGoingAverage / iteration_, 2);

        double standard_error = sqrt(variance / iteration_);

        standard_error_ = standard_error;

        return onGoingAverage / iteration_;
    }
    double PriceWithoutVarience()
    {
        double price = 0;
        for (int current_cycle = 0; current_cycle < iteration_; current_cycle++)
        {
            price += OneCycle();
        }

        return price / iteration_;
    }
};

#endif