#include "localCoordinateChange.h"
#include "constants.h"
#include "VectorField.h"
#include "readSections.h"

fromSectionCoordinates::fromSectionCoordinates(IVector q_,IMatrix A_)
{
	H = IMap(EnergyFormula());
	H.setParameter("mu",mu);
	K = IMap(Energy2BPFormula());
	K.setParameter("mu",mu);

	g = IMap("var:u,s,alpha,I,y4,y5,eps;fun:u,s,alpha,y4,y5,0,eps;");

	A=A_;
	q=q_;
	h=H(q0())[0];
	kappa0=K(q0())[0];
}

//   G(w1,w2,alpha,I,y4,y5,eps) =
//      (H(q+A(w1,w2,alpha,y4,y5,0,eps))-h,
//		 K(q+A(w1,w2,alpha,y4,y5,0,eps))-I-kappa0)
//
//   G(z) =
//      (H(q+A*g(z))-h,
//		 K(q+A*g(z))-I-kappa0).
IVector fromSectionCoordinates::G(const IVector &z) const
{
	IVector r(2);

	// z=(w1,w2,alpha,I,y4,y5,eps)
	interval I=z[3];
	
	r[0] = H(q+A*g(z))[0]-h;
	r[1] = K(q+A*g(z))[0]-I-kappa0;
	return r;
}

IVector Z(const IVector &x,const IVector &y)
{
	IVector z(7);
	// x=(w1,w2,alpha,I,eps)
	// y=(y4,y5)
	// z=(w1,w2,alpha,I,y4,y5,eps)
	for(int i=0;i<4;i++) z[i]=x[i];
	z[4] = y[0];
	z[5] = y[1];
	z[6] = x[4];
	return z;
}

IVector fromSectionCoordinates::G(const IVector &x,const IVector &y) const
{
	return G(Z(x,y));
}

IMatrix fromSectionCoordinates::DG(const IVector &z) const
{
	IMatrix dG(2,7);
	
	IMatrix dH=H[q+A*g(z)]*A*g[z];
	IMatrix dK=K[q+A*g(z)]*A*g[z];
	
	// G(z) =
	//      (H(q+A*g(z))-h,
	//		 K(q+A*g(z))-I-kappa0).
	// z=(w1,w2,alpha,I,y4,y5,eps)
	
	for(int i=0;i<7;i++) dG[0][i] = dH[0][i];
	for(int i=0;i<7;i++) dG[1][i] = dK[0][i];
	dG[1][3] = dK[0][3] - 1.0; // derivative with respect to I.
	return dG;
}

IMatrix fromSectionCoordinates::dGdx(const IVector &z) const
{
	IMatrix dG=DG(z);
	// x=(w1,w2,alpha,I,eps)
	// z=(w1,w2,alpha,I,y4,y5,eps)
	IMatrix D(2,5);
	for(int i=0;i<2;i++)
		for(int j=0;j<4;j++)
			D[i][j] = dG[i][j];
	D[0][4] = dG[0][6];
	D[1][4] = dG[1][6];
	return D;
}

IMatrix fromSectionCoordinates::dGdy(const IVector &z) const
{
	IMatrix dG=DG(z);
	// x=(w1,w2,alpha,I,eps)
	// z=(w1,w2,alpha,I,y4,y5,eps)
	IMatrix D(2,2);
	
	D[0][0] = dG[0][4]; D[0][1] = dG[0][5];
	D[1][0] = dG[1][4]; D[1][1] = dG[1][5];
	return D;
}

IMatrix fromSectionCoordinates::dGdy(const IVector &x,const IVector &y) const
{
	return dGdy(Z(x,y));
}

IMatrix fromSectionCoordinates::dGdx(const IVector &x,const IVector &y) const
{
	return dGdx(Z(x,y));
}

IVector fromSectionCoordinates::Y(const IVector &x) const
{
	IVector y(2);
	IVector x0=midVector(x);
	// computing an initial approximation:
	for(int i=0;i<20;i++) y = midVector(y - gauss(dGdy(x0,y),G(x0,y)));
	
	// enlarging the initial approximation:
	for(int i=0;i<3;i++)
	{
		y = midVector(y) - gauss(dGdy(x,y),G(x,midVector(y)));
	}
	IVector y0=midVector(y);
	y = midVector(y) + 1.5*(y-y0);

	// validation using Interval Newton method:
	IVector N = y0 - gauss(dGdy(x,y),G(x,y0));
	if(subsetInterior(N,y)) return N;
	cout << "computation of the local change of coordinates failed. Aborting." << endl;
	abort();
	return IVector(2);
}

IMatrix fromSectionCoordinates::DY(const IVector &x) const
{
	IVector y=Y(x);
	return -gaussInverseMatrix(dGdy(x,y))*dGdx(x,y);
}

IVector fromSectionCoordinates::omega(const IVector &x) const
{
	return g(Z(x,Y(x)));
}

IMatrix fromSectionCoordinates::Domega(const IVector &x) const
{
	IMatrix dy=DY(x);
	IMatrix D(7,5);
	// x=(w1,w2,alpha,I,eps)
	// p=g(x,Y(x))=(w1,w2,alpha,y4,y5,0,eps)
	for(int i=0;i<3;i++) D[i][i]=1.0;
	for(int i=0;i<2;i++)
		for(int j=0;j<5;j++)
			D[i+3][j] = dy[i][j];
	D[6][4]=1.0;
	return D;
}

///////////////////

fromLocalCoordinates::fromLocalCoordinates(IVector q,IMatrix A,IMatrix B_)
{
	psi=fromSectionCoordinates(q,A);
	B=B_;
}

///////////////////

toSectionCoordinates::toSectionCoordinates(IVector q_,IMatrix A)
{
	q=q_;
	Ainv=gaussInverseMatrix(A);
	K = IMap(Energy2BPFormula());
	K.setParameter("mu",mu);
	kappa0=K(q0())[0];
}

////////////////////////////////
// The change of coordinates computes 
// 		x = (w1,w2,alpha,I,eps). 
// From a point on the surface of section. This is done as follows. The
//    Ainv(p-q) = (w1,w2,alpha,y4,y5,0,eps)
// so w1,w2,alpha and eps simply follow from projections from Ainv(p-q).
// The coordinate I is K(p)-kappa0, so taking I=K(p)-kappa0
IVector toSectionCoordinates::image(IVector p) const
{
	IVector x(5);
	p=Ainv*(p-q);
	for(int i=0;i<3;i++) x[i]=p[i]; // w1,w2,alpha
	x[3]=K(p)[0]-kappa0; // I
	x[4]=p[6]; // eps
	return x;
}

/////////////////////////////
// Looking at the comments for above toSectionCoordinates::image() 
// we can see that the derivative is the composition of 
//   - Ainv
//   - derivative of the projections
//   - derivative of K.
IMatrix toSectionCoordinates::derivative(const IVector &p) const
{
	IMatrix C(5,7);
	C[0][0]=1; // w1
	C[1][1]=1; // w2
	C[2][2]=1; // alpha
	C[4][6]=1; // eps
	IMatrix D=C*Ainv;
	IMatrix dK=K[p];
	for(int i=0;i<7;i++) D[3][i]=dK[0][i]; 
	return D;
}

////////////////////////////////

toLocalCoordinates::toLocalCoordinates(IVector q,IMatrix A,IMatrix B)
{
	psiInv=toSectionCoordinates(q,A);
	Binv=gaussInverseMatrix(B);
}


