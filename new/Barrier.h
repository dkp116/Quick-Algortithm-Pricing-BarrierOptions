#ifndef BARRIER_H
#define BARRIER_H

#include "Option.h"

class DownAndOut : public Option {
private:
    double barrier_;

public:
    DownAndOut(ExerciseType exerciseType, OptionType optionType, double strike, double barrier)
        : Option(strike, exerciseType, optionType), barrier_(barrier) {}
        
    double payoff(double currentValue) const override;
};

#endif // BARRIER_H
