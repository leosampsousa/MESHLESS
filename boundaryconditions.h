#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include"node.h"
#include<vector>
#include"sparsedefinition.h"
using namespace std;


class BoundaryConditions
{
public:
    BoundaryConditions();
    virtual~BoundaryConditions();
    virtual void impose (StiffnessSparseMatrix* stiffness, ForceSparseVector* force, const vector<Node>& constrainedNodes, int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f,bool *ok) = 0;
};

#endif // BOUNDARYCONDITIONS_H
