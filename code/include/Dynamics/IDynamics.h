#ifndef IDYNAMICS_H
#define IDYNAMICS_H


class IDynamics{
    
    public:
    virtual ~IDynamics() = default;
    virtual double evolve(double TimeIncrement) =0 ;
};


#endif