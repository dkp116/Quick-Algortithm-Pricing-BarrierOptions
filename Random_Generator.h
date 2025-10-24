#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>
#include <thread>
#include <functional>

class RandomGenerator{
    private:
    static unsigned seed(){
        static std::random_device rd;
        return rd()^ static_cast<unsigned>(std::hash<std::thread::id>{}(std::this_thread::get_id()));
    }
    public:
    static std::mt19937& getGenerator(){
        thread_local static std::mt19937 gen(seed());
        return gen;
    }
};

#endif