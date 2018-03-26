#ifndef CONSTRAINEDDOF_H
#define CONSTRAINEDDOF_H
#include "dof.h"
class constrainedDOF : public DOF{
private:
    double restriction;
public:
    constrainedDOF(int equationNumber,double restriction): DOF(equationNumber),restriction(restriction){}

    bool isConstrained(){return true;}

    double getRestriction() const {return restriction;}
};






#endif // CONSTRAINEDDOF_H
