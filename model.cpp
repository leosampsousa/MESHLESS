#include "model.h"
#include "constraineddof.h"

StiffnessSparseMatrix Model::getStiffness() const
{
    return stiffness;
}

ForceSparseVector Model::getForce() const
{
    return force;
}

WeightFunction *Model::getWeightFunction() const
{
    return weightFunction;
}

const vector<Node>& Model::getNodes() const
{
    return nodes;
}

Model::Model(int n_node, double min, double max, double alpha, double beta):
    min(min),max(max)
{
#if PRINT
    cout << "Criando Modelo"<<endl;
#endif
    this->weightFunction=new Spline();
    double step = (max-min)/(n_node-1);
    DOF *dof;
    for(int i = 0;i<n_node;i++){
        if (i == 0 || i == n_node-1){
            dof = new constrainedDOF (i,0);
        }else{
            dof = new DOF(i);
        }
        nodes.push_back(Node(i+1, step*i, alpha*step ,beta*step,dof));
    }
    stiffness.resize(n_node,n_node);
    stiffness.setZero();
    force.resize(n_node,1);
    force.setZero();

}

Model::~Model()
{
    nodes.clear();
}

void Model::calculateSystem(MethodId method, int polynomialOrder, const Load &T, Load &K, Load &f, int pInt,int n_node)
{
    //erro era que um vetor de <CoefficientTriplet estava sendo passado como argumento da função calculateLocalSystem que pedia
    //apenas um CoefficientTriplet
    //vou mudar a definição da função
    vector<CoefficientTriplet>coeffK, coeff;
#if PRINT
    cout<<"\t\tcalculando sistema"<<endl;
#endif

    for(Node node : nodes)
    {
        node.calculateLocalSystem(this->min,this->max,method,polynomialOrder,this->weightFunction,nodes,T,K,f,&coeffK,&coeff, pInt,n_node);
    }

    stiffness.setFromTriplets(coeffK.begin(),coeffK.end());
    force.setFromTriplets(coeff.begin(),coeff.end());
}

VectorXd Model::update(VectorXd u_hat, int polynomialOrder)
{
#if PRINT
    cout<< "u chapeu atualizado"<<endl;
#endif

    VectorXd result(u_hat.rows());
    for(Node node : nodes){
        result(node.getDof()->getEquationNumber()) = node.calculateDeslocamento(u_hat, polynomialOrder,this->getWeightFunction(),nodes);
    }

    return result;
}


vector<Node> Model::getConstrainedNodes()
{
#if PRINT
    cout<<"pegando nós com restrição"<<endl;
#endif
    vector<Node> constrainedNodes;
    for (Node node : nodes){
        if(node.getDof()->isConstrained()){
            constrainedNodes.push_back(node);
        }
    }

    return constrainedNodes;
}
