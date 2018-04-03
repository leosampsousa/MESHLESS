#include "analyzer.h"
#include <iomanip>

VectorXd Analyzer::getU() const
{
    return u;
}

VectorXd Analyzer::getU_hat() const
{
    return u_hat;
}

Analyzer::Analyzer(int polynomialOrder, Model model, Solver solver, MethodId method, BoundaryConditions* boundaryconditions, Load T, Load K, Load f)
    : polynomialOrder(polynomialOrder), model(model), solver(solver), method(method), boundaryconditions(boundaryconditions), T(T), K(K), f(f)
{
#if PRINT
    cout<<"criando analyzer"<<endl;
#endif
}

Analyzer::~Analyzer()
{
    delete boundaryconditions;
}

void Analyzer::analyze(int pInt, int n_node)
{
#if PRINT
    cout<<"analisando o modelo"<<endl;
#endif
    this->model.calculateSystem(this->method,this->polynomialOrder,this->T,this->K,this->f,pInt, n_node);
#if PRINT
    cout<<"pós calculo sistema"<<endl;
#endif
    //função model tem q receber um ponteiro para stifness


   // MatrixXd stiffness = model.getK();

    StiffnessSparseMatrix stiffness = model.getStiffness();
#if PRINT
    cout <<"recebendo a Matriz ANTES\n"<<MatrixXd(stiffness)<<endl;
#endif
    //VectorXd force = model.getF();

    ForceSparseVector force = model.getForce();
#if PRINT
    cout <<"recebendo o vetor ANTES\n"<<MatrixXd(force)<<endl;
#endif
    this->boundaryconditions->impose(&stiffness,&force,model.getConstrainedNodes(),polynomialOrder,model.getWeightFunction(),model.getNodes(),this->T,this->K,this->f);

#if PRINT
    cout <<"\n\nrecebendo a Matriz DEPOIS\n"<<MatrixXd(stiffness)<<endl;
    cout <<"recebendo o vetor DEPOIS\n"<<MatrixXd(force)<<endl;

    cout<<"impondo condições de contorno"<<endl;
#endif
    u_hat = this->solver.solve(stiffness,force);

#if PRINT
    cout<<"resolvendo sistema"<<endl;
    cout<<"u chapeu: \n"<< scientific << setprecision(15) <<u_hat<<endl;
#endif
    //adicionar a cada nó seu deslocamento
    u = model.update(u_hat,this->polynomialOrder);
    cout<<"u: \n"<< scientific << setprecision(15) <<u<<endl;
#if PRINT
    cout<<"atualizando modelo"<<endl;
#endif
}


