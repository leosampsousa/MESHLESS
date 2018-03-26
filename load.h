#ifndef LOAD_H
#define LOAD_H

#include <iostream>
using namespace std;
class Load
{
private:
    double constant = 1;

public:
    Load(double constant);
    ~Load();
    double calculateContribuition() const;
};

#endif // LOAD_H
