#include "Barrier.h"
#include <algorithm>

double DownAndOut::payoff(double currentValue) const {
    if (currentValue <= barrier_) return 0.0;
    if(optionType_ == OptionType::Call)
        return std::max(currentValue - strike_, 0.0);
    else
        return std::max(strike_ - currentValue, 0.0);
}
