#if USE_INTERFACE

#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}

#else

#include "model.h"
#include "solver.h"
#include "methodid.h"
#include "boundaryconditions.h"
#include"analyzer.h"
#include "lagrangemultiplier.h"
#include "direct.h"
#include "solucaoexata.h"
int main()
{
    int pInt = 5;
    int polynomialOrder=2; //esta com a diferença de 1
    int n_node=50;
    double min=0;
    double max=3;
   // int pInt=2;
    double alpha = 0.9;
    double beta= 2.1;
    double passoAlpha = 0.1;
    double passoBeta = 0.1;
    VectorXd u;
    VectorXd solucaoExata(n_node);
    MethodId method = MethodId::DMLPG1;
   // Model model(n_node, min, max, alpha, beta);
    Solver solver;

    SolucaoExata solucao;

    Load T(20);
    Load K(2);
    Load f(10);
/*
    Model model(n_node, min, max, alpha, beta);
    BoundaryConditions* boundaryConditions = new Direct();
    Analyzer analyzer(polynomialOrder,model,solver, method,boundaryConditions,T,K,f);
    analyzer.analyze(pInt,n_node);
*/

    int contador = 0;
    for(double i = 0; i<=max ; i+=max/(n_node-1)){
        solucaoExata(contador) = solucao.calculate(i);
        contador++;
    }


    VectorXd melhor_u;
    double melhorAlpha;
    double melhorBeta;
    int melhorpInt;
    double error;
    double comparador = 100;
    for(pInt = 3; pInt<=6;pInt++){
        for(alpha = 0.6;alpha <=2.4;alpha+=passoAlpha){
            for(beta = 1.2; beta <=3.0;beta+=passoBeta){
                if(alpha>beta) continue;
                VectorXd norma;
                Model model(n_node, min, max, alpha, beta);
                BoundaryConditions* boundaryConditions = new Direct();
                Analyzer analyzer(polynomialOrder,model,solver, method,boundaryConditions,T,K,f);
                analyzer.analyze(pInt,n_node);
                u = analyzer.getU();
                norma = solucaoExata - u;
                error = (norma.norm())/(solucaoExata.norm());
                if(error < comparador){
                    comparador = error;
                    melhorAlpha = alpha;
                    melhorBeta = beta;
                    melhorpInt = pInt;
                    melhor_u = u;
                    //cout<<"O menor erro foi : "<<comparador<<"\nele foi obtido por meio do alpha igual a: "<<melhorAlpha<<" e do beta igual a: "<<melhorBeta<<endl;
                }
                cout <<"\n-------------"<< endl;
                cout <<"pInt = "<<pInt << endl;
                cout <<"alpha = " <<alpha << endl;
                cout <<"beta = " <<beta << endl;
                cout <<"erro = " << error << endl;
            }
        }
    }


    cout<<"\nO menor erro foi : "<<comparador<<"\nele foi obtido por meio do alpha igual a: "<<melhorAlpha<<" e do beta igual a: "<<melhorBeta<<" e com "<< melhorpInt <<" pontos de integração"<<endl;
    cout <<"u: \n"<<melhor_u <<endl;
    cout <<"solucaoExata: \n"<<solucaoExata <<endl;

    return 0;
}


#endif
