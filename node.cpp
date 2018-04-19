#include "node.h"
#include "integration.h"
#include "monomialbases.h"
#include "constraineddof.h"
#include <iomanip>
#include <fstream>


DOF *Node::getDof() const
{
    return dof;
}

void Node::setDof(DOF *value)
{
    dof = value;
}

int Node::getId() const
{
    return id;
}

double Node::getDeslocamento() const
{
    return deslocamento;
}

double Node::calculateDeslocamentoQualquer(double position, VectorXd u_hat, int polynomialOrder, WeightFunction* weight, const vector<Node> &nodes, bool *ok)
{
    double deslocamento;
    int cont=0;
#if PRINT
        cout<<"\t\t\tCalculando deslocamento final"<<endl;
#endif

    DMLPG_parameters* parameters = calculateShapeFunctionParameters(position,polynomialOrder,weight,nodes, ok);
    if(!*ok) return 0.0;
   // DMLPG_parameters* parameters = calculateShapeFunctionParameters(3/2,polynomialOrder,weight,nodes);
    for(Node node:parameters->adjacentNodes){
        deslocamento += node.calculateShapeFunction(parameters, cont++)*u_hat(node.dof->getEquationNumber());
    }
    delete parameters;

    return deslocamento;
}

void Node::shapeGraphic(int totalPoints ,int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, bool *ok, ofstream *file, int n_node, int min , int max , MethodId method,Load T, Load K, Load f)
{
    DMLPG_parameters *parameters;
    int N = 2*totalPoints;
    double minLocal, maxLocal,step;
    int count;
    double shape[n_node];
    double shapeDiff[n_node];
    double sumShape[n_node];


    *file << "\"node_"<<id<<"_point"<<"\" \"node_"<<id<<"_accumulator"<<"\" " ;
    *file << "\"node_"<<id<<"_accumulatorDiff"<<"\" " << "\"node_"<<id<<"_sum"<<"\" ";

    for(int i = 0 ; i<n_node;i++){

        *file << "\"node_"<<id<<"_shape_"<<i<<"\" \"node_"<<id<<"_shapeDiff_"<<i<<"\" " ;
        *file << "\"node_"<<id<<"_sumShape_"<<i<<"\" " ;
    }

    *file<<endl;
    if(isLeft() || isRight(n_node)){
        N = totalPoints;
    }

    minLocal = position - weightRadius - 0.1;
    maxLocal = position + weightRadius + 0.1;

    if(minLocal < min){
        minLocal = min;
    }

    if(maxLocal > max){
        maxLocal = max;
    }
    step = (maxLocal - minLocal)/N;


    for(int i = 0 ; i <= N ; i++){
        count = 0;
        double point = minLocal + i*step;

        for(int i = 0;i < n_node;i++){
            shape[i]=0;
            shapeDiff[i]=0;
            sumShape[i]=0;
        }

        parameters = calculateIntegrationPointsParameters(point,method,polynomialOrder,weight,nodes,T,K,f,ok);
        if(!*ok){
            return;
        }
        double accumulator=0;
        double accumulatorDiff = 0;
        double sum =0;
        for(Node node:parameters->adjacentNodes){
            double shapeValue = (parameters->testFunction*node.calculateShapeFunction(parameters, count))*parameters->K;
            double shapeDiffValue = (parameters->testeFunctionDiff*node.calculateShapeFunctionDerivative(parameters, method, count))*parameters->T;
            shape[node.id-1]= shapeValue;
            shapeDiff[node.id-1]= shapeDiffValue;
            sumShape[node.id-1] = shapeValue + shapeDiffValue;
            accumulator += shapeValue;
            accumulatorDiff += shapeDiffValue;
            sum += shapeValue + shapeDiffValue;
            count++;
        }

        *file << point << " " <<accumulator<<" "<< accumulatorDiff << " " << sum;

        for(int i = 0 ; i<n_node;i++){

           *file << " " << shape[i] << " " <<shapeDiff[i]<< " " <<sumShape[i];
        }

        *file<<endl;


        delete parameters;

    }


}

double Node::calculateDeslocamento(VectorXd u_hat, int polynomialOrder, WeightFunction* weight, const vector<Node> &nodes,bool *ok)
{
    deslocamento = calculateDeslocamentoQualquer(position,u_hat,polynomialOrder,weight,nodes,ok);
    if(!*ok){
        return 0.0;
    }
    return deslocamento;
}

Node::Node(int id, double position,  double testRadius, double weightRadius, DOF *dof) :
    id(id), position(position), testRadius(testRadius), weightRadius(weightRadius), deslocamento(0), dof(dof)
{
#if PRINT
    cout<<"\tCriando nó " << id<<" "<<position<<endl;
#endif
}

Node::~Node()
{

}

DMLPG_parameters *Node :: calculateShapeFunctionParameters(double point,int polynomialOrder, WeightFunction *weight,const vector<Node> &nodes, bool *ok){

    DMLPG_parameters* parameters = new DMLPG_parameters();
    MonomialBases monomialBases;
    parameters->p = monomialBases.calculate(point,polynomialOrder);

    VectorXd w(nodes.size());

    int adjacentNodes_size=0;
    for (const Node node :  nodes){
       double valor = weight->calculate(point,node.position,node.weightRadius);
       if (valor > 0){
           w(adjacentNodes_size++)=valor;
           parameters->adjacentNodes.push_back(node);
       }
    }
    w.conservativeResize(adjacentNodes_size);

    MatrixXd P (adjacentNodes_size,polynomialOrder);
    for(int i=0;i<adjacentNodes_size;i++){
        P.row(i)= monomialBases.calculate(parameters->adjacentNodes[i].position,polynomialOrder).transpose();
    }

    MatrixXd B = P.transpose()*w.asDiagonal()*w.asDiagonal();


    MatrixXd A = B*P;


    if(A.fullPivLu().isInvertible()){
        *ok = true;
#if PRINT
        cout<<"\t\t\t\té invertivel"<<endl;
#endif
    }else{
        *ok = false;
        cout<<"\t\t\t\tnão é invertivel"<<endl;
        return NULL;

    }
    MatrixXd A1 = A.inverse();

    parameters->A1B = A1*B;
    return parameters;
}

DMLPG_parameters *Node::calculateIntegrationPointsParameters(double point, MethodId method, int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f,bool *ok)
{
    DMLPG_parameters* parameters;
    if (method == MethodId::DMLPG1){
        parameters = new DMLPG_parameters();
    }else{
        parameters = new MLPG_parameters();
    }
    MonomialBases monomialBases;

    parameters->T = T.calculateContribuition();
    parameters->K = K.calculateContribuition();
    parameters->f = f.calculateContribuition();
    parameters->testFunction = weight->calculate(point,position,testRadius);
    parameters->testeFunctionDiff = weight->differentiate(point,position,testRadius);
    parameters->p = monomialBases.calculate(point,polynomialOrder);
    parameters->pDiff = monomialBases.differentiate(point,polynomialOrder);


 #if PRINT
    cout<<"\n\n\t\t\t\tcalculando os parametros do ponto de integração " << point<<" e posição: "<<position<<endl;

    cout<<"\t\t\t\tspline pro ponto: "<< scientific << setprecision(15)<< parameters->testFunction<<endl;
    cout<<"\t\t\t\tderivada spline pro ponto: "<< parameters->testeFunctionDiff<<endl;

    cout<<"\t\t\t\tbase monomial ponto: "<< parameters->p.transpose()<<endl;
    cout<<"\t\t\t\tderivada base monomial ponto: "<< parameters->pDiff.transpose()<<endl;
#endif

     VectorXd w(nodes.size());

    int adjacentNodes_size=0;
    for (const Node node :  nodes){
        double valor = weight->calculate(point,node.position,node.weightRadius);
        if (valor > 0){
            w(adjacentNodes_size++)=valor;
            parameters->adjacentNodes.push_back(node);
        }
    }
    w.conservativeResize(adjacentNodes_size);

    MatrixXd P (adjacentNodes_size,polynomialOrder);
    for(int i=0;i<adjacentNodes_size;i++){
        P.row(i)= monomialBases.calculate(parameters->adjacentNodes[i].position,polynomialOrder).transpose();
    }

#if PRINT
    cout<<"\t\t\t\tw :"<<MatrixXd (w.asDiagonal())<<endl;
    cout<<"\t\t\t\tP: "<<P<<endl;
#endif
   // MatrixXd B = P.transpose()*w.asDiagonal()*w.asDiagonal();
    MatrixXd B = P.transpose()*w.asDiagonal();


    MatrixXd A = B*P;


    if(A.fullPivLu().isInvertible()){
        *ok = true;
#if PRINT
        cout<<"\t\t\t\té invertivel"<<endl;
#endif
    }else{
        *ok = false;
        cout<<"\t\t\t\tnão é invertivel"<<endl;
        return NULL;

    }
    MatrixXd A1 = A.inverse();

    parameters->A1B = A1*B;
#if PRINT
    cout<<"\t\t\t\tmatriz A: "<<A<<endl;
    cout<<"\t\t\t\tmatriz A-1: "<<A.inverse()<<endl;
    cout<<"\t\t\t\tmatriz B: "<<B<<"\n\n"<<endl;

#endif
     if (method == MethodId::MLPG1){
         //utilizAR O W BARRA
         VectorXd w_barraDif(adjacentNodes_size);

         for(int i=0;i<adjacentNodes_size;i++){
            //usar spline.diferentiate
            double valor = weight->differentiate(point,parameters->adjacentNodes[i].position,weightRadius);
            w_barraDif(i)=valor;
           // w_barraDif(i)=valor*valor;
         }

         //derivada do A, Derivada do B
         MatrixXd B_diff = P.transpose()*w_barraDif.asDiagonal();
         MatrixXd A_diff = B_diff*P;

         //formulas 36 e 37
         //depois da soma desconsiderando o p transposto guardar no A1diff

        //nao consegui adicionar a parameters
         ((MLPG_parameters*)(parameters))->A1BDiff = (-A1*A_diff*A1)*B + A1*B_diff;

     }

    //calcula w e w' e add em adjacentNodes
    return parameters;

}

void Node::calculateContribuitionK(DMLPG_parameters *parameters, vector<CoefficientTriplet> *coeffK,int j, Node node, double w, MethodId method)
{

    //o I é do nó corrente , o J é do no parametro, e o valor é eq 66 e 67
    //equation number vem do dof
    double contribuitionK;
    double equationNumberI, equationNumberJ;
    equationNumberI = dof->getEquationNumber();
    equationNumberJ = node.getDof()->getEquationNumber();
    contribuitionK =
            w*(((parameters->testeFunctionDiff)*(node.calculateShapeFunctionDerivative(parameters,method,j))*parameters->T) +
            (parameters->K*parameters->testFunction*node.calculateShapeFunction(parameters,j)));

#if PRINT
    cout<<"\t\t\t\t\tshape function do nó: "<<node.id<<" é: "<<node.calculateShapeFunction(parameters,j)<<endl;
    cout<<"\t\t\t\t\tderivada da shape function do nó: "<<node.id<<" é: "<<node.calculateShapeFunctionDerivative(parameters,method,j)<<endl;

    cout<<"\t\t\t\t\tcontribuiação K("<<equationNumberI<<","<<equationNumberJ<<"): "<<contribuitionK<<endl;
#endif
    coeffK->push_back(CoefficientTriplet(equationNumberI, equationNumberJ, contribuitionK));
    //cout<<"\t\t\t\t\tcontribuiação K("<<equationNumberI<<","<<equationNumberJ<<"): "<<coeffK<<endl;

}

void Node::calculateContribuitionf(DMLPG_parameters *parameters, vector<CoefficientTriplet> *coefff, double w)
{
    double contribuitionf;
    double equationNumberI;
    equationNumberI = dof->getEquationNumber();
    contribuitionf =
            w*(parameters->f*parameters->testFunction);
#if PRINT
    cout<<"\t\t\t\tvalores de f: "<<contribuitionf<<"\n\n"<<endl;
#endif

    coefff->push_back(CoefficientTriplet(equationNumberI, 0, contribuitionf));
}

void Node::calculateContribuitionLeft(DMLPG_parameters *parameters, vector<CoefficientTriplet> *coeffK, int j, Node node, MethodId method)
{
    double contribuitionLeft;
    double equationNumberI, equationNumberJ;
    equationNumberI = dof->getEquationNumber();
    equationNumberJ = node.getDof()->getEquationNumber();

    contribuitionLeft = (-parameters->T)*(node.calculateShapeFunctionDerivative(parameters,method,j));
    coeffK->push_back(CoefficientTriplet(equationNumberI, equationNumberJ, contribuitionLeft));
}

void Node::calculateContribuitionRight(DMLPG_parameters *parameters, vector<CoefficientTriplet> *coeffK, int j, Node node, MethodId method)
{
    double contribuitionRight;
    double equationNumberI, equationNumberJ;
    equationNumberI = dof->getEquationNumber();
    equationNumberJ = node.getDof()->getEquationNumber();

    contribuitionRight = (parameters->T)*(node.calculateShapeFunctionDerivative(parameters,method,j));
    coeffK->push_back(CoefficientTriplet(equationNumberI, equationNumberJ, contribuitionRight));
}

void Node::calculateLocalSystem(double min, double max, MethodId method, int polynomialOrder,
                                WeightFunction *weight, const vector<Node> &nodes,
                                const Load &T,  Load &K,  Load &f,
                                vector<CoefficientTriplet> *coeffK, vector<CoefficientTriplet> *coefff, int pInt,int n_node,bool *ok)
{
#if PRINT
    cout<<"\n----------------------------\n\t\tcalculando sistema local do nó: "<<id<<endl;
#endif

    Integration integration;
    for(int i=0;i<pInt;i++){
        double intPosition = 0;
        double intWeight = 0;
        integration.createIntegrationPoint(this->id,i,pInt,min,max,testRadius, &intWeight, &intPosition);


        if (isLeft()){
            intPosition = (testRadius/2.0)*(1 + intPosition);
            intWeight = -(testRadius/2.0)*intWeight;
        }else if(isRight(n_node)){
            intPosition = ((max - min) - testRadius/2) + (testRadius/2)*intPosition;
            intWeight = - (testRadius/2.0)*intWeight;
        }else{
            intPosition = this->position + testRadius*intPosition;
            intWeight = -intWeight*testRadius;
        }
#if PRINT
        cout<<"\t\t\tPosição ponto de integração "<<i << ": "<<intPosition<<endl;
        cout<<"\t\t\tPeso ponto de integração "<<i << ": "<<intWeight<<endl;
#endif
        DMLPG_parameters* parameters = calculateIntegrationPointsParameters(intPosition,method,polynomialOrder,weight,nodes,T,K,f,ok);
        if(!*ok){
            delete parameters;
            return;
        }

        int cont=0;
        for(Node node:parameters->adjacentNodes){
#if PRINT
            cout<<"\t\t\t\tcalculando contribuição K do nó adjacente: "<<node.id<<endl;
#endif
            calculateContribuitionK(parameters,coeffK,cont++,node,intWeight,method);
        }
        calculateContribuitionf(parameters,coefff,intWeight);


        delete parameters;
    }
#if PRINT
    cout<<"\t\t\tDepois dos pontos de integração"<<"\n\n"<<endl;
#endif

    //se nó esquedo
    int cont=0;
    if (isLeft()){
#if PRINT
        cout<<"\t\t\tRecalculando parametros para o nó esquerdo"<<endl;
#endif
        DMLPG_parameters* parameters = calculateIntegrationPointsParameters(min,method,polynomialOrder,weight,nodes,T,K,f,ok);
        if(!*ok){
            delete parameters;
            return;
        }
        for(Node node:parameters->adjacentNodes){
            calculateContribuitionLeft(parameters,coeffK,cont++,node,method);
        }
        delete parameters;
    }



    //se nó direito
    if(isRight(n_node)){
#if PRINT
        cout<<"\t\t\tRecalculando parametros para o nó direito"<<endl;
#endif
        DMLPG_parameters* parameters = calculateIntegrationPointsParameters(max,method,polynomialOrder,weight,nodes,T,K,f,ok);
        cont=0;
        if(!*ok){
            delete parameters;
            return;
        }
        for(Node node:parameters->adjacentNodes){
            calculateContribuitionRight(parameters,coeffK,cont++,node,method);
        }
        delete parameters;
    }



}



double Node::calculateShapeFunction(DMLPG_parameters *parameters, int j)
{

    return (parameters->p.transpose()*parameters->A1B)(0,j);
}

double Node::calculateShapeFunctionDerivative(DMLPG_parameters *parameters, MethodId method, int j)
{
    //MatrixXd shapeFunctionDerivative;
    if (method == MethodId::DMLPG1){
        //cout <<parameters->pDiff<<endl;
       // cout << parameters->A1B<<endl;

        return (parameters->pDiff.transpose()*parameters->A1B)(0,j);
    }else {
        return (parameters->pDiff.transpose()*parameters->A1B + parameters->p.transpose()*((MLPG_parameters*)(parameters))->A1BDiff)(0,j);
    }

}

bool Node::isLeft()
{
   if(this->id==1){
       return true;
   }
   return false;
}

bool Node::isRight(int n_node)
{
    if(this->id==n_node){
        return true;
    }
    return false;
}

double Node::calculateImposeContribution(DMLPG_parameters *parameters, int j, Node node){
    double contribuition = node.calculateShapeFunction(parameters,j);
    return contribuition;
}

void Node :: imposeCondition(StiffnessSparseMatrix *stiffness, ForceSparseVector *force, int polynomialOrder, WeightFunction *weight, const vector<Node> &nodes, const Load &T, const Load &K, const Load &f, bool *ok){
    DMLPG_parameters* parameters = calculateIntegrationPointsParameters(position, MethodId::DMLPG1,polynomialOrder, weight, nodes, T, K, f, ok);
    if(!*ok){
        delete parameters;
        return;
    }
#if PRINT
    cout<<"posição: "<<position<<endl;
#endif
    int cont=0;
    double valor;
    for(Node node:parameters->adjacentNodes){
        valor = calculateImposeContribution(parameters,cont++,node);
#if PRINT
        cout<<"novo valor de stifiness depois de imposicao: "<<valor<<endl;
#endif
        stiffness->coeffRef(dof->getEquationNumber(),node.dof->getEquationNumber()) = valor;
    }

    force->coeffRef(dof->getEquationNumber(),0) = ((constrainedDOF*)dof)->getRestriction();
    delete parameters;

}


