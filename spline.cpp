#include "spline.h"

Spline::Spline()
{
#if PRINT
    cout<<"\tCriando spline"<<endl;
#endif
}

Spline::~Spline()
{

}
// x_i é o valor onde ela é centrada !
double Spline::calculate(double x, double x_i, double radius)
{
    double d_i;
    d_i = fabs(x - x_i);
    if(d_i >= radius){
        return 0;
    }
    return 1 - 6 * pow(d_i/radius,2) + 8* pow(d_i/radius,3) - 3* pow(d_i/radius,4);
}

double Spline::differentiate(double x, double x_i, double radius)
{
    double d_i;
    double p_i;
    d_i = fabs(x - x_i);
    p_i = x_i-x;

    return ((12*pow(d_i,2))/pow(radius,4))*p_i - ((24*d_i)/pow(radius,3))*p_i + (12/pow(radius,2))*p_i;

}
