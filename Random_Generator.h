#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>

class RandomGenerator{
    public:
    static std::mt19937& getGenerator(){
        static std::mt19937 gen(std::random_device{}());
        return gen;
    }
};

#endif