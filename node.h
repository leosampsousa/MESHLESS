#ifndef NODE_H
#define NODE_H
#include"weightfunction.h"
#include"methodid.h"
#include"load.h"
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include"spline.h"
#include"dof.h"
#include"sparsedefinition.h"

using namespace Eigen;
using namespace std;
class Node;

struct DMLPG_parameters {
    double K;
    double T;
    double f;
    double testFunction;
    double testeFunctionDiff;
    VectorXd p;
    MatrixXd A1B;
    VectorXd pDiff;
    vector<Node> adjacentNodes;
};

struct MLPG_parameters : public DMLPG_parameters {
    MatrixXd A1BDiff;
};

typedef Triplet<double> CoefficientTriplet;
class Node
{
private :
    int id;
    double position;
    double testRadius;
    double weightRadius;
    double deslocamento;
    DOF *dof;

public:
    Node(int id,double position, double testRadius, double weightRadius,DOF *dof);
    ~Node();
    void calculateLocalSystem(double min, double max, MethodId method, int polynomialOrder, WeightFunction* weight, const vector<Node>& nodes, const Load &T, Load &K, Load &f, vector<CoefficientTriplet>* coeffK, vector<CoefficientTriplet>* coefff, int pInt, int n_node);
    DMLPG_parameters *calculateIntegrationPointsParameters(double point, MethodId method, int polynomialOrder, WeightFunction* weight, const vector<Node>& nodes, const Load &T, const Load &K, const Load &f);
    //faltando o "parameters"
    static DMLPG_parameters* calculateShapeFunctionParameters(double point,int polynomialOrder, WeightFunction* weight,const vector<Node> &nodes);
    void calculateContribuitionK(DMLPG_parameters* parameters, vector<CoefficientTriplet>* coeffK, int j, Node node, double w, MethodId method);
    void calculateContribuitionf(DMLPG_parameters* parameters, vector<CoefficientTriplet>* coefff, double w);
    void calculateContribuitionLeft(DMLPG_parameters* parameters, vector<CoefficientTriplet>* coeffK, int j, Node node, MethodId method);
    void calculateContribuitionRight(DMLPG_parameters* parameters, vector<CoefficientTriplet>* coeffK, int j, Node node, MethodId method);
    double calculateShapeFunction(DMLPG_parameters* parameters, int j);
    double calculateShapeFunctionDerivative(DMLPG_parameters* parameters, MethodId method, int j);
    bool isLeft();
    bool isRight(int n_node);
    void imposeCondition(StiffnessSparseMatrix *stiffness, ForceSparseVector *force, int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f);
    double calculateImposeContribution(DMLPG_parameters *parameters, int j, Node node);
    DOF *getDof() const;
    void setDof(DOF *value);
    int getId() const;
    double getDeslocamento() const;
    double calculateDeslocamento(VectorXd u_hat, int polynomialOrder, WeightFunction* weight, const vector<Node> &nodes);
    static double calculateDeslocamentoQualquer(double position, VectorXd u_hat, int polynomialOrder, WeightFunction* weight, const vector<Node> &nodes);
};



#endif // NODE_H
