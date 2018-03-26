#ifndef MONOMIALBASES_H
#define MONOMIALBASES_H
#include<Eigen/Core>
#include<vector>
#include <iostream>
using namespace std;
using namespace Eigen;

class MonomialBases
{
public:
    MonomialBases();
    ~MonomialBases();
    VectorXd calculate(double x,int polynomialOrder);
    //se for 2d, adicionar um vector<VectorXd>
    VectorXd differentiate (double x, int polynomialOrder);
};

#endif // MONOMIALBASES_H
