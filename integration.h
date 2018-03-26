#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <iostream>
#include <vector>
using namespace std;
class Integration
{
public:
    Integration();
    ~Integration();
    vector<int> id;
    void createIntegrationPoint(int id, int i , int numXi, double x0, double xN, double radius, double *_w, double *_p );
};

#endif // INTEGRATION_H
