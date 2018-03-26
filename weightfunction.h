#ifndef WEIGHTFUNCTION_H
#define WEIGHTFUNCTION_H

class WeightFunction{
public:
    virtual double calculate(double x,double x_f, double radius)=0;
    virtual double differentiate(double x, double x_f, double radius)=0;
};

#endif // WEIGHTFUNCTION_H
