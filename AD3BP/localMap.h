#ifndef localMap_h
#define localMap_h

#include "capd/capdlib.h"
#include <iostream>
using namespace capd;
using namespace std;

#include "localCoordinateChange.h"

////////////////////////////////////
// Class: localMap
// 
// Purpose: This class allows us to consider two surfaces of sections
// whose local coordinates are based on a choice of respectively:
//    q0, A0, B0
// and 
//    q1, A1, B1.
// See "5.2 Choices of local coordinates on sections" from the paper,
// and also documentation from localCoordinateChange.h/cpp.
// For a point 
//    x0 = (u0,s0,alpha0,I0,eps0)
// on the first surface of section, the local map passes this point
// to the coordinates alpha,J,r2,phi2,R2,PHI2,eps, computes a Poincare map P
// to the next section, and then passes to the local coordinates 
// u1,s1,alpha1,I1,eps1 on that section.
//
// The resulting map is:
//     f = psiInv o P o psi
// 
// Members:
//     psi    - change from coordinates u0,s0,alpha0,I0,eps0 to alpha,J,r2,phi2,R2,PHI2,eps
//              on the first section given by q0,A0,B0.
//     psiInv - change of coordinates from section1 in the original coordinates alpha,J,r2,phi2,R2,PHI2,eps
//              to the local coordinates (u1,s1,alpha1,I1,eps1) on section1.
//     solver - this is a technical item. It is used by CAPD to initiate a Poincare map to section1.
//     section1 - this is the surface of section which is orthogonal to the fifth column of A1 at q1.
//     P      - map which takes a point from the state space and integrates it to section1.
//     K      - the Kepler part of the Hamiltonian in coordinates alpha,J,r2,phi2,R2,PHI2,eps.
//              The K is used by us to compute the Keplerian part of the energy, which is I+kappa0.
//     q      - two points q[0],q[1] at which the surfaces of section are positioned.
//     A      - two matrices A[0], A[1]. The fifthe column of A[1] determines section1.
//              The matrices are used for the computation of local coordinates on the sections.
//     B      - two matrices B[0], B[1] which are used for the computation of the local coordinates
//              on the surfaces of section.
//
// The only finction of interest is:
//     image() - this is the function of main interest. It takes a point
//               X0 = (u0,s0,alpha0,I0,eps0)
//               and computes
//               X1 = (u1,s1,alpha1,I1,eps1)
//               by passing to the coordinates alpha,J,r2,phi2,R2,PHI2,eps, flowing to section1
//               and then passing to the local coordinates on section1.
//               To speed up the computations in other parts of the program 
//               this function computes several things:
//                  fx0 - we name the midpoint of the box as X0 and call it x0. 
//                        the fx0 will be a bound on the image of x0.
//                  Df  - This will be the bound on the derivative of f on the whole initial set X0.
//                  dI  - This is the bound on the change change of
//                        		K(P(x))-K(x)
//                        This is computed using (64) from the paper and enxures that
//                              K(P(x))-K(x) \in eps*dI.          
//               The function returns X1=f(X0), and passes fx0, Df, dI by reference.
class localMap
{
private:
	fromLocalCoordinates psi;
	toLocalCoordinates psiInv;
	
	IOdeSolver* solver;
	IAffineSection* section1;
	IPoincareMap* P;
	IMap* K;

	vector<IVector> q;
	vector<IMatrix> A;
	vector<IMatrix> B;
public:	
	localMap(IVector q0,IVector q1,IMatrix A0,IMatrix A1,IMatrix B0,IMatrix B1,IOdeSolver &solver_);
	
	localMap(){}
	~localMap();
	localMap& operator=(const localMap& other);
	
	IVector image(const IVector &X0,IVector &fx0,IMatrix &Df,interval &dI) const;
	IVector operator()(const IVector &X0,IVector &fx0,IMatrix &Df,interval &dI) const {return image(X0,fx0,Df,dI);}
	
	vector<IVector> get_q() const {return q;}
	vector<IMatrix> get_A() const {return A;}
	vector<IMatrix> get_B() const {return B;}
};

vector<localMap> sequenceOfLocalMaps(IOdeSolver &solver);

#endif