#include "direct.h"

Direct::Direct()
{

}

Direct::~Direct()
{

}

void Direct::impose(StiffnessSparseMatrix *stiffness, ForceSparseVector *force, const vector<Node> &constrainedNodes,int polynomialOrder,WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f)
{
    for(Node node : constrainedNodes){
        node.imposeCondition(stiffness,force,polynomialOrder,weight,nodes,T,K,f);
    }
}
