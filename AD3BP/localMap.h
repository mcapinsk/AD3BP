#ifndef localMap_h
#define localMap_h

#include "capd/capdlib.h"
#include <iostream>
using namespace capd;
using namespace std;

#include "localCoordinateChange.h"

class localMap
{
private:
	interval h;
	
	fromLocalCoordinates psi;
	toLocalCoordinates psiInv;
	
	IMap* f;
	IOdeSolver* solver;
	IAffineSection* section;
	IPoincareMap* P;
	IMap* K;

public:	

	vector<IVector> q;
	vector<IMatrix> A;
	vector<IMatrix> B;
	
	localMap(IVector q0,IVector q1,IMatrix A0,IMatrix A1,IMatrix B0,IMatrix B1,interval h,IMap &f_,IOdeSolver &solver_);
	localMap(){}
	~localMap();
	localMap& operator=(const localMap& other);
	
	IVector image(const IVector &x,IVector &fx0,IMatrix &Df,interval &dI) const;
	IVector operator()(const IVector &x,IVector &fx0,IMatrix &Df,interval &dI) const {return image(x,fx0,Df,dI);}
	
	vector<IVector> get_q() const {return q;}
	vector<IMatrix> get_A() const {return A;}
	vector<IMatrix> get_B() const {return B;}
	interval get_h() const {return h;}
};

vector<localMap> sequenceOfLocalMaps(IMap &f,IOdeSolver &solver);

#endif