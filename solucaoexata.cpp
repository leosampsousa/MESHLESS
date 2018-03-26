#include "solucaoexata.h"
#include <math.h>

SolucaoExata::SolucaoExata()
{

}

double SolucaoExata::calculate(double valor)
{
    return (5 * pow(M_E,(-valor/sqrt(10))) * (pow(M_E, 3/sqrt(10)) - pow(M_E, valor/sqrt(10))) * (pow(M_E, valor/sqrt(10)) - 1))
            /
            (1 + pow(M_E, 3/sqrt(10)));
}
