#ifndef localCoordinateChange_h
#define localCoordinateChange_h

#include "capd/capdlib.h"
using namespace capd;
using namespace std;
using namespace capd::matrixAlgorithms;

//////////////////////////////////
// Class: fromSectionCoordinates
//
// Purpose: The aim of this class is to pass from local coordinates
// to a section attached at a point q, and spanned by vectors given by
// the matrix A. (See "5.2 Choices of local coordinates on sections" 
// from the paper.) We take a point x from the local coordinates and compute
// its image in the full state space of the system:
//     q + A*omega(x).
// 
//
// As discussed in "5.2 Choices of local coordinates on sections" 
// from the paper, the local coordinates on the surface of section are:
//   w1,w2,alpha,I,eps
// they are mapped to the surface of section as
//   x = q+A(w1,w2,alpha,y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps),0,eps)
// where y4,y5 are functions for which
//   H(q+A(w1,w2,alpha,y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps),0,eps)) = h.
//   K(q+A(w1,w2,alpha,y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps),0,eps)) = K(q0)+I = kappa0+I.
//
// For a given (w1,w2,alpha,I,eps) y4,y5 are found by solving
//   H(q+A(w1,w2,alpha,y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps),0,eps)) - h = 0
//   K(q+A(w1,w2,alpha,y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps),0,eps)) - kappa0 - I = 0
// using Newton's method.
//
// Note: The zero appearing in the formulae is associated with the fact that the 
// coordinates are chosen on a Poincare section which is of the form
//   S = {q+A(w1,w2,w3,w4,w5,0,eps)}
// The vector
//   v = A*(0,0,0,0,0,1,0)
// is chosen to be orthogonal to 
//   A*ei  
// for all i different than 6, where      i
// ei = (0,..,0,1,0,..,0)
//
// Remark: the fact that v = A*(0,0,0,0,0,1,0) is orthogonal to 
// all other vectors is rigorously ensured by 
//    makeOrthogonalToFifthColumn()
// which is called in getLinearChanges().
// 
// We consider an auxilary function: G:R^7 -> R^2
//   G(w1,w2,alpha,I,y4,y5,eps) =
//      (H(q+A(w1,w2,alpha,y4,y5,0,eps))-h,
//		 K(q+A(w2,w2,alpha,y4,y5,0,eps))-I-kappa0.
//
// We will also consider an auxilary function: g:R^7 -> R^7
//    g(w1,w2,alpha,I,y4,y5,eps) = (w1,w2,alpha,y4,y5,0,eps)
// which allows us to write
//   G(y) = G(w1,w2,alpha,I,y4,y5,eps) =
//      (H(q+A*g(y))-h,
//		 K(q+A*g(y))-I-kappa0.
//
// DG will be a function returning a 2 x 7 matrix, the derivative of G.
//
// We will also have a function dGdx returning a 2 x 5 matrix
// that is the derivative of G with respect to (w1,w2,alpha,I,eps).
//
// We will also have a function dGdy returning a 2 x 2 matrix
// that is the derivative of G with respect to (y4,y5).
//
// We have also a function Y:R^5 -> R^2 defined as
//    Y(w1,w2,alpha,I,eps) = (y4(w1,w2,alpha,I,eps),y5(w1,w2,alpha,I,eps))=(y4,y5)
// where (y4,y5) are obtained by solving G=0.
//
// We define 
//   omega(x) = omega(w1,w2,alpha,I,eps)
//      = g(w1,w2,alpha,Y[0](x),Y[1](x),eps)
//      = (w1,w2,alpha,Y[0](x),Y[1](x),0,eps)
// 
// With this notation our coordinate change can be written as
//    q + A*omega(x)
class fromSectionCoordinates
{
private:
	IMatrix A;
	IVector q;
	IMap H,K; // H is the Hamiltonian of the full 3bp. 
			  // K is the Kepler part of the Hamiltonian.

	// g is an auxiliary function		   
	IMap g;   // g(w1,w2,alpha,I,y4,y5,eps) = (w1,w2,alpha,y4,y5,0,eps)
	interval h; // this is H(q0)
	interval kappa0; // this is K(q0)
	
public:
	
	// our convention for the notation used for the variables is that:
	//   x will be 5 dimensional x=(w1,w2,alpha,I,eps)
	//   y will be 2 dimensional y=(y4,y5) 
	//   z will be 7 dimensional z=(w1,w2,alpha,I,y4,y5,eps)

	IVector G(const IVector &z) const;
	IVector G(const IVector &x,const IVector &y) const;

	IMatrix DG(const IVector &z) const;
	IMatrix dGdx(const IVector &z) const;
	IMatrix dGdy(const IVector &z) const;
	IMatrix dGdy(const IVector &x,const IVector &y) const;
	IMatrix dGdx(const IVector &x,const IVector &y) const;

	IVector Y(const IVector &x) const;
	IMatrix DY(const IVector &x) const;

	IVector omega(const IVector &x) const;
	IMatrix Domega(const IVector &x) const;
	
	IVector image(const IVector &x) const {return q+A*omega(x);}
	IMatrix derivative(const IVector &x) const {return A*Domega(x);}
	
	fromSectionCoordinates(IVector q,IMatrix A);
	fromSectionCoordinates(){}
	
	IVector operator()(const IVector &x) const {return image(x);}

	// Note: our convention that the [] operator returns the derivative.
	// This is the standing convention used throughout CAPD.
	IMatrix operator[](const IVector &x) const {return derivative(x);}

	// In CAPD the integrator takes the so called C0Rect2Set
	// as an initial value for the integration. The C0Rect2Set
	// is a set represented as 
	//     q+A*R
	// where 
	//     q - is a point (IVector)
	//     A - is a matrix (IMatrix)
	//     R - is a box (IVector)
	// The constructor for a C0Rect2Set takes (q,A,R).
	C0Rect2Set C0Set(const IVector &x) const 
	{
		IVector B=omega(x);
		IVector B0=midVector(B);
		// q+A*omega(x) = q+A*B = q+A*B0+A*(B-B0) 
		return C0Rect2Set(q+A*B0,A,B-B0);
	}
	// The C1Rect2Set plays the same role as the C0Rect2Set
	// and has the same representation. The difference between the
	// two is that C1Rect2Set is passed to a solver if C1 computations
	// are required.
	C1Rect2Set C1Set(const IVector &x) const 
	{
		IVector B=omega(x);
		IVector B0=midVector(B);
		return C1Rect2Set(q+A*B0,A,B-B0);
	}
};

///////////////////////////////
// Class: fromLocalCoordinates
//
// This class equips the class fromSectionCoordinates with
// an additional linear change of coordinates, given by a matrix B.
//
// The local coordinates on the surface of section are
//    (w1,w2,alpha,I,eps) = B*(u,s,alpha,I,eps).
//
// The total change of coordinates is 
//     psi(B*x)
// where 
//     x = (u,s,alpha,I,eps)
// and 
//     psi is an object of class fromSectionCoordinates.
class fromLocalCoordinates
{
private:
	fromSectionCoordinates psi;
	IMatrix B;
public:
	IVector image(const IVector &x) const {return psi(B*x);}
	C0Rect2Set C0Set(const IVector &x) const {return psi.C0Set(B*x);}
	C1Rect2Set C1Set(const IVector &x) const {return psi.C1Set(B*x);}

	// Note our convention that psi[p] returns the derivative of 
	// psi at p. In other words, psi[p] = Dpsi(p).
	// This is the standing convention used throughout CAPD.
	IMatrix derivative(const IVector &x) const {return psi[B*x]*B;}
	
	IVector operator()(const IVector &x) const {return image(x);}
	IMatrix operator[](const IVector &x) const {return derivative(x);}
	
	fromLocalCoordinates(IVector q,IMatrix A,IMatrix B);
	fromLocalCoordinates(){}
};

///////////////////////////////////////
// Class: toSectionCoordinates
// 
// This class gives the inverse change to the fromSectionCoordinates.
// In other words, it passes from the state space coordinates 
//    p = q + A*omega(x)
// to the local coordinates 
//    x = (w1,w2,alpha,I,eps). 
// This is done as follows. The
//    Ainv(p-q) = (w1,w2,alpha,y4,y5,0,eps)
// so w1,w2,alpha and eps simply follow from projections.
// The coordinate I is K(p)-kappa0, so taking I=K(p)-kappa0 allows us to compute x.
class toSectionCoordinates
{
private:
	IMatrix Ainv;
	IVector q;
	IMap K;
	interval kappa0;
public:
	toSectionCoordinates(IVector q,IMatrix A);
	toSectionCoordinates(){}
	IVector image(IVector p) const;
	IMatrix derivative(const IVector &p) const;
	
	IVector operator()(const IVector &p) const {return image(p);}
	IMatrix operator[](const IVector &p) const {return derivative(p);}
};

class toLocalCoordinates
{
private:
	toSectionCoordinates psiInv;
	IMatrix Binv;
public:
	IVector image(const IVector &x) const {return Binv*psiInv(x);}
	IMatrix derivative(const IVector &x) const {return Binv*psiInv[x];}
	
	IVector operator()(const IVector &x) const {return image(x);}
	IMatrix operator[](const IVector &x) const {return derivative(x);}
	
	toLocalCoordinates(IVector q,IMatrix A,IMatrix B);
	toLocalCoordinates(){}
};

#endif