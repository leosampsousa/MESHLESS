#include "solver.h"


Solver::Solver()
{
#if PRINT
    cout<<"chamando solver"<<endl;
#endif
}

Solver::~Solver()
{

}

VectorXd Solver::solve(StiffnessSparseMatrix stiffness,ForceSparseVector force)
{
    VectorXd u_hat;
    SparseLU<StiffnessSparseMatrix>lu;
    stiffness.makeCompressed();
    lu.analyzePattern(stiffness);
    lu.factorize(stiffness);
    u_hat = lu.solve(force);

    return u_hat;

}
