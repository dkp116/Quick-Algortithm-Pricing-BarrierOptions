#ifndef OPTION_H
#define OPTION_H

enum class ExerciseType {American,European};
enum class OptionType {Call, Put};
class Option{

    public:
    virtual ~Option() = default;
    virtual double Payoff(double currentValue) const = 0;
    double GetRebate(){return rebate_;}
    double GetStrike(){return strike_;}

    protected:
    Option(ExerciseType exerciseType , OptionType optionType, double strike , double rebate) : 
        strike_(strike) , exerciseType_(exerciseType) , optionType_(optionType), rebate_(rebate) {}
        
    double strike_;
    ExerciseType exerciseType_;
    OptionType optionType_;
    double rebate_;     //Should we have the price of a vanilla call here? 
};


#endif