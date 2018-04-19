#ifndef ANALYZER_H
#define ANALYZER_H
#include"model.h"
#include"solver.h"
#include"methodid.h"
#include"boundaryconditions.h"

class Analyzer
{
private:
    int polynomialOrder;
    Model model;
    Solver solver;
    MethodId method;
    BoundaryConditions* boundaryconditions;
    Load T;
    Load K;
    Load f;
    VectorXd u;
    VectorXd u_hat;
public:
    Analyzer(int polynomialOrder, Model model, Solver solver, MethodId method, BoundaryConditions* boundaryconditions, Load T, Load K, Load f);
    ~Analyzer();
    bool analyze(int pInt, int n_node);
    VectorXd getU() const;
    VectorXd getU_teste() const;
    VectorXd getU_hat() const;
};

#endif // ANALYZER_H
