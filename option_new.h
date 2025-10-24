#ifndef OPTIONNEW_H
#define OPTIONNEW_H

enum class ExerciseType {American,European};
enum class OptionType {call, put};
class Option{

    public:
    virtual ~Option() = default;

    protected:
    Option(double strike , ExerciseType exerciseType , OptionType optionType) : 
        strike_(strike) , exerciseType_(exerciseType) , optionType_(optionType) {}
        
    double strike_;
    ExerciseType exerciseType_;
    OptionType optionType_;

};




#endif