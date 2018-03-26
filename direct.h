#ifndef DIRECT_H
#define DIRECT_H
#include"boundaryconditions.h"

class Direct:public BoundaryConditions
{
public:
    Direct();
    ~Direct();
    void impose (StiffnessSparseMatrix* stiffness, ForceSparseVector* force, const vector<Node>& constrainedNodes, int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f) override;
};

#endif // DIRECT_H
