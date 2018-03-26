#ifndef SPLINE_H
#define SPLINE_H
#include "weightfunction.h"
#include <iostream>
#include<math.h>
using namespace std;
class Spline : public WeightFunction
{
public:
    Spline();
    ~Spline();
    double calculate(double x, double x_i, double radius) override;
    double differentiate(double x, double x_i, double radius) override;
};

#endif // SPLINE_H
