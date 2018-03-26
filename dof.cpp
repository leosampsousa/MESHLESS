#include"dof.h"

DOF::DOF(int equationNumber):equationNumber(equationNumber){}

double DOF::getEquationNumber()
{
    return equationNumber;
}

bool DOF::isConstrained()
{
    return false;
}
