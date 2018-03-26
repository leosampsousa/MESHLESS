#ifndef DOF_H
#define DOF_H
class DOF{
private:
    int equationNumber;
public:
    DOF(int equationNumber);
    double getEquationNumber();
    virtual bool isConstrained();
};

#endif // DOF_H
