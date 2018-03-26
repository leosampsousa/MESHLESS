#ifndef SOLVER_H
#define SOLVER_H
#include<Eigen/Core>
#include <iostream>
#include"sparsedefinition.h"
using namespace Eigen;
using namespace std;
class Solver
{
public:
    Solver();
    ~Solver();
    VectorXd solve(StiffnessSparseMatrix stiffness,ForceSparseVector force);
};

#endif // SOLVER_H
