#ifndef localCoordinateChange_h
#define localCoordinateChange_h

#include "capd/capdlib.h"
using namespace capd;
using namespace std;
using namespace capd::matrixAlgorithms;

// local coordinates are:
//   u,s,alpha,I,eps
// they are mapped to the section as
//   x = q+A(u,s,alpha,v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps),0,eps)
// where v0,v1 are functions for which
//   H(q+A(u,s,alpha,v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps),0,eps)) = h.
//   K(q+A(u,s,alpha,v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps),0,eps)) = K(q)+I.
//
// For a given (u,s,alpha,I,eps) v0,v1 are found by solving
//   H(q+A(u,s,alpha,v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps),0,eps)) - h = 0
//   K(q+A(u,s,alpha,v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps),0,eps)) - K(q) - I = 0
// using Newton's method.
//
// Note: The zero appearing in the formulae is associated with the fact that the 
// coordinates are chosen on a Poincare section which is of the form
//   S = {q+A(x1,x2,x3,x4,x5,0,eps)}
// The vector
//   v = A*(0,0,0,0,0,1,0)
// is chosen to be orthogonal to 
//   A*ei  
// for all i different than 6, where 
//              i
// ei = (0,..,0,1,0,..,0)
// 
// We consider an auxilary function: G:R^7 -> R^2
//   G(u,s,alpha,I,v0,v1,eps) =
//      (H(q+A(u,s,alpha,v0,v1,0,eps))-h,
//		 K(q+A(u,s,alpha,v0,v1,0,eps))-I-K(q)).
//
// We will consider an auxilary function: g:R^7 -> R^7
//    g(u,s,alpha,I,v0,v1,eps) = (u,s,alpha,v0,v1,0,eps)
// which allows us to write
//   G(y) = G(u,s,alpha,I,v0,v1,eps) =
//      (H(q+A*g(y))-h,
//		 K(q+A*g(y))-I-K(q)).
// DG will be a function returning a 2 x 7 matrix, the derivative of G.
// We will have a function dGdx returning a 2 x 5 matrix
// that is the derivative of G with respect to (u,s,alpha,I,eps).
// We will also have a function dGdv returning a 2 x 2 matrix
// that is the derivative of G with respect to (v0,v1).
//
// We have also a function V:R^5 -> R^2 defined as
//    V(u,s,alpha,I,eps) = (v0(u,s,alpha,I,eps),v1(u,s,alpha,I,eps))=(v0,v1)
// where (v0,v1) are obtained by solving G=0.
//
// We define 
//   omega(x) = omega(u,s,alpha,I,eps)
//      = (u,s,alpha,V[0](x),V[1](x),0,eps)
// 
// With this notation our coordinate change can be written as
//    q + A*omega(x)
class fromSectionCoordinates
{
private:
	IMatrix A;
	IVector q;
	IMap H,K; // added to the class so that the constructor is not called 
			  // each time that V is computed
	IMap g;   // g(u,s,alpha,I,v0,v1,eps) = (u,s,alpha,v0,v1,0,eps)
	interval h; // this is H(q)
	interval Kq; // this is K(q)

	IMap to_R;
	
public:
	interval I1(const IVector &x) const;
	IMatrix dI1(const IVector &x) const;
	
	// auxhiliary functions for I2:
	// here f(u,s,alpha,I,c1,eps)=H(q+A(u,s,alpha,I,c1,0,eps))-h
	// our convention is that: 
	//   x will be 5 dimensional x=(u,s,alpha,I,eps)
	//   y will be 6 dimensional y=(u,s,alpha,I,c1,eps)
	//   p will be 7 dimensional p=(u,s,alpha,I,c1,0,eps)
	// our convention is that:
	//   v will be 2 dimensional v=(v0,v1) 
	//   x will be 5 dimensional x=(u,s,alpha,I,eps)
	//   y will be 7 dimensional y=(u,s,alpha,I,v0,v1,eps)

	IVector G(const IVector &y) const;
	IVector G(const IVector &x,const IVector &v) const;

	IMatrix DG(const IVector &y) const;
	IMatrix dGdx(const IVector &y) const;
	IMatrix dGdv(const IVector &y) const;
	IMatrix dGdv(const IVector &x,const IVector &v) const;
	IMatrix dGdx(const IVector &x,const IVector &v) const;

	IVector V(const IVector &x) const;
	IMatrix DV(const IVector &x) const;

	IVector omega(const IVector &x) const;
	IMatrix Domega(const IVector &x) const;
	
	IVector image(const IVector &x) const {return q+A*omega(x);}
	IMatrix derivative(const IVector &x) const {return A*Domega(x);}
	
	fromSectionCoordinates(IVector q,IMatrix A,interval h);
	fromSectionCoordinates(){}
	
	IVector operator()(const IVector &x) const {return image(x);}
	IMatrix operator[](const IVector &x) const {return derivative(x);}

	C0Rect2Set C0Set(const IVector &x) const 
	{
		IVector B=omega(x);
		IVector B0=midVector(B);
		return C0Rect2Set(q+A*B0,A,B-B0);
	}
	C1Rect2Set C1Set(const IVector &x) const 
	{
		IVector B=omega(x);
		IVector B0=midVector(B);
		return C1Rect2Set(q+A*B0,A,B-B0);
	}
};

class fromLocalCoordinates
{
private:
	fromSectionCoordinates psi;
	IMatrix B;
public:
	IVector image(const IVector &x) const {return psi(B*x);}
	C0Rect2Set C0Set(const IVector &x) const {return psi.C0Set(B*x);}
	C1Rect2Set C1Set(const IVector &x) const {return psi.C1Set(B*x);}

	IMatrix derivative(const IVector &x) const {return psi[B*x]*B;}
	
	IVector operator()(const IVector &x) const {return image(x);}
	IMatrix operator[](const IVector &x) const {return derivative(x);}
	
	fromLocalCoordinates(IVector q,IMatrix A,interval h);
	fromLocalCoordinates(IVector q,IMatrix A,IMatrix B,interval h);
	fromLocalCoordinates(){}
};

// class toSectionCoordinates
// {
// private:
// 	IMatrix Ainv;
// 	IVector q;
// public:
// 	toSectionCoordinates(IVector q,IMatrix A);
// 	toSectionCoordinates(){}
// 	IVector image(IVector p) const;
// 	IMatrix derivative(const IVector &p) const;
	
// 	IVector operator()(const IVector &p) const {return image(p);}
// 	IMatrix operator[](const IVector &p) const {return derivative(p);}
// };

class toSectionCoordinates
{
private:
	IMatrix Ainv;
	IVector q;
	IMap K;
	interval Kq;
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
	
	toLocalCoordinates(IVector q,IMatrix A);
	toLocalCoordinates(IVector q,IMatrix A,IMatrix B);
	toLocalCoordinates(){}
};

#endif