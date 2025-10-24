#ifndef BARRIER_H
#define BARRIER_H

#include "option_new.h"

class DownAndOut : public Option {
private:
    double barrier_;

public:
    DownAndOut(double strike, double barrier, ExerciseType exerciseType, OptionType optionType)
        : Option(strike, exerciseType, optionType), barrier_(barrier) {}
        
    double payoff(double currentValue) const override;
};

#endif // BARRIER_H
