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

#define LOOP_MAIN false

#include "model.h"
#include "solver.h"
#include "methodid.h"
#include "boundaryconditions.h"
#include"analyzer.h"
#include "lagrangemultiplier.h"
#include "direct.h"
#include "solucaoexata.h"
#include <fstream>
#include <stdio.h>
#include <sstream>
int main()
{
    //teste
    ofstream file("output_MLPG_5.txt");
    ofstream fileShape("output_MLPG_GRAPHIC_3.txt");
    int pInt = 3;
    int polynomialOrder=2; //esta com a diferença de 1
    int n_node=5;
    double min=0;
    double max=3;
    double alpha = 0.55;
    double beta= 1.63;
    double passoAlpha = 0.01;
    double passoBeta = 0.01;
    double deslocamentoPontoCentral;
    VectorXd u;
    VectorXd solucaoExata(n_node);
    MethodId method = MethodId::MLPG1;
   // Model model(n_node, min, max, leoalpha, beta);
    Solver solver;

    SolucaoExata solucao;

    Load T(20);
    Load K(2);
    Load f(10);

#if !LOOP_MAIN

    Model model(n_node, min, max, alpha, beta);
    bool ok = true;
    model.shapeGraphics(30,polynomialOrder,&ok,&fileShape,n_node,method,T,K,f);
    if(!ok){
        cout<<"error"<<endl;
    }

    BoundaryConditions* boundaryConditions = new Direct();
    Analyzer analyzer(polynomialOrder,model,solver, method,boundaryConditions,T,K,f);
    analyzer.analyze(pInt,n_node);


#else

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
    for(pInt = 3; pInt<=5;pInt++){
        for(alpha = 0.4;alpha <=1.0001;alpha+=passoAlpha){

            file << "\"Beta_"<< "n" << n_node << "_" << "m" << polynomialOrder<<"_"<<"i"<<pInt <<"_" <<"a"<< alpha << "\" ";
            file << "\"Error_"<< "n" << n_node << "_" << "m" << polynomialOrder<<"_"<<"i"<<pInt <<"_" <<"a"<< alpha << "\" ";
            file << "\"u_central_"<< "n" << n_node << "_" << "m" << polynomialOrder<<"_"<<"i"<<pInt <<"_" <<"a"<< alpha << "\""<<endl;


            for(beta = 0.8; beta <=2.0001;beta+=passoBeta){
                if(alpha>beta) continue;
                VectorXd norma;
                Model model(n_node, min, max, alpha, beta);
                BoundaryConditions* boundaryConditions = new Direct();
                Analyzer analyzer(polynomialOrder,model,solver, method,boundaryConditions,T,K,f);
                bool ok = analyzer.analyze(pInt,n_node);
                if(!ok) continue;
                deslocamentoPontoCentral = Node::calculateDeslocamentoQualquer(max/2.0,analyzer.getU_hat(),polynomialOrder,model.getWeightFunction(),model.getNodes(),&ok);
                if(!ok){
                    continue;
                }
                u = analyzer.getU();
                norma = solucaoExata - u;
                error = (norma.norm())/(solucaoExata.norm());

                file << beta<< " "<<error << " "<< deslocamentoPontoCentral<<endl;

                if(error < comparador){
                    comparador = error;
                    melhorAlpha = alpha;
                    melhorBeta = beta;
                    melhorpInt = pInt;
                    melhor_u = u;
                    //cout<<"O menor erro foi : "<<comparador<<"\nele foi obtido por meio do alpha igual a: "<<melhorAlpha<<" e do beta igual a: "<<melhorBeta<<endl;
                }
                /*
                cout <<"\n-------------"<< endl;
                cout <<"pInt = "<<pInt << endl;
                cout <<"alpha = " <<alpha << endl;
                cout <<"beta = " <<beta << endl;
                cout <<"erro = " << error << endl;
                */
            }
        }
    }



    file.close();

    cout<<"\nO menor erro foi : "<<comparador<<"\nele foi obtido por meio do alpha igual a: "<<melhorAlpha<<" e do beta igual a: "<<melhorBeta<<" e com "<< melhorpInt <<" pontos de integração"<<endl;
    cout <<"u: \n"<<melhor_u <<endl;
    cout <<"solucaoExata: \n"<<solucaoExata <<endl;





#endif
    return 0;
}
#endif
