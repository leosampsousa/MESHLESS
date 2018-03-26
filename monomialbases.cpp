#include "monomialbases.h"

MonomialBases::MonomialBases()
{
#if PRINT
    cout<<"\t\t\tcriando base monomial"<<endl;
#endif
}

MonomialBases::~MonomialBases()
{
}

VectorXd MonomialBases::calculate(double x, int polynomialOrder){

    VectorXd vector(polynomialOrder);
    for(int i = 0; i<polynomialOrder;i++){
        vector(i)=pow(x,i);
    }
    return vector;

}

VectorXd MonomialBases::differentiate(double x, int polynomialOrder){
    VectorXd vector(polynomialOrder);
    vector(0)=0;
    for(int i = 0; i<polynomialOrder-1;i++){
        vector(i + 1)=(i+1)*pow(x,i);
    }
    return vector;
}
