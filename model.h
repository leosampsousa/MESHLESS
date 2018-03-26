#ifndef MODEL_H
#define MODEL_H
#include"methodid.h"
#include"load.h"
#include<vector>
#include"node.h"
#include<math.h>
using namespace std;
#include"sparsedefinition.h"


class Model
{
private:
    StiffnessSparseMatrix stiffness;
    ForceSparseVector force;
    double min;
    double max;
    void createIntegrationPoints();
    vector<Node> nodes;
    WeightFunction* weightFunction;
public:
    Model(int n_node, double min, double max,double alpha,double beta);
    ~Model();
    void calculateSystem(MethodId method, int polynomialOrder, const Load &T, Load &K, Load &f, int pInt, int n_node);
    VectorXd update(VectorXd hat_u, int polynomialOrder);
    //MatrixXd getK();
   //  VectorXd getF();
    vector<Node> getConstrainedNodes();

    StiffnessSparseMatrix getStiffness() const;
    ForceSparseVector getForce() const;
    WeightFunction *getWeightFunction() const;
   const vector<Node> &getNodes() const;
};

#endif // MODEL_H
