#include "load.h"

Load::Load(double constant):constant(constant)
{
#if PRINT
    cout<<"criando load"<<endl;
#endif
}

Load::~Load()
{

}

double Load::calculateContribuition() const
{
#if PRINT
    cout<<"\t\t\tcalculando contribuicao de load"<<endl;
#endif
    return constant;
}
