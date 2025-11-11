#ifndef OPTION_H
#define OPTION_H

enum class ExerciseType {American,European};
enum class OptionType {Call, Put};
class Option{

    public:
    virtual ~Option() = default;
    virtual double payoff(double currentValue) const = 0;

    protected:
    Option(double strike , ExerciseType exerciseType , OptionType optionType) : 
        strike_(strike) , exerciseType_(exerciseType) , optionType_(optionType) {}
        
    double strike_;
    ExerciseType exerciseType_;
    OptionType optionType_;

};


#endif