#ifndef BARRIER_H
#define BARRIER_H

#include "Option.h"

class DownAndOut : public Option {
private:
    double barrier_;

public:
    DownAndOut(ExerciseType exerciseType, OptionType optionType, double strike, double barrier, double rebate)
        : Option(exerciseType, optionType,strike, rebate), barrier_(barrier) {}

    double GetBarrier(){return barrier_;}

    double  Payoff(double currentValue) const override;


};

#endif // BARRIER_H
