#ifndef readOrbits_h
#define readOrbits_h

#include "capd/capdlib.h"
#include <iostream>
using namespace capd;
using namespace std;

vector<IVector> readOrbit();
vector<IMatrix> getLinearChanges();
vector<IMatrix> getLocalLinearChanges();

#endif