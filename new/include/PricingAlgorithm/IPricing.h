#ifndef IPricing_h
#define IPricing_h

#include "Dynamics/IDynamics.h"
#include "Stock/Stock.h"
#include "Options/Option.h"
#include <iostream>
#include <chrono>

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

    template <typename Func>
    double ProfileUsingTime(Func func)
    {
        auto start = std::chrono::high_resolution_clock::now();

        double result = func();

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        time_ = elapsed.count();

        return result;
    }

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
        auto pricing = [&]()
        {
            if (includeStandardError_ == StandardErrorCalculation::NotIncluded)
                return PriceWithoutVariance();
            else
                return PriceWithVariance();
        };

        if (includeTime_ == Time::Included)
            return ProfileUsingTime(pricing);

        return pricing();
    }

    double PriceWithVariance()
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
    double PriceWithoutVariance()
    {
        double price = 0;
        for (int current_cycle = 0; current_cycle < iteration_; current_cycle++)
        {
            price += OneCycle();
        }

        return price / iteration_;
    }

    double GetTime(){
        return time_;
    }

    double GetStandardError(){
        return standard_error_;
    }
};

#endif