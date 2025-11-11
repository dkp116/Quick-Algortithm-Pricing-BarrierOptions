#ifndef BROWNIANBRIDGE_H
#define BROWNIANBRIDGE_H
#include "IPricing.h"
#include "Stock.h"

class BrownianBridge : public IPricing{
    private:
    double iteration_;
    public:
    BrownianBridge(std::shared_ptr<Stock> stock , std::shared_ptr<Option> option , double iteration) : IPricing(stock, option) , iteration_(iteration) {}
    double Price() override;

};
//change the constuctor
//make the shared pointer to only that, but this defeats the point of the factor method right 
//only see that as a vaible solution here:
//make a shared ptr of the specific solution here and set them equal , then use the "local shared ptr to calcualte"
//what if we want 2 solutions ??? 
// not sure how to over come this, unless when we call price(option_type) and then assign the ptr?
//shared_ptr and just assign it in consturction, but how do we know which one to assign too?
//identifier so when called it gives the dynamcis and this is used to dertmine which one we will use like a swithc statement right? 
//can we have multiple constuctors as this does not need to be generic right and i am pretty sure that 

//take in specifc child ptrs and have the interface take in the adult? can i still pass the child in??? google when you have internet

#endif