#include "localCoordinateChange.h"
#include "constants.h"
#include "VectorField.h"

fromSectionCoordinates::fromSectionCoordinates(IVector q_,IMatrix A_,interval h_)
{
	H = IMap(EnergyFormula());
	H.setParameter("mu",mu);
	K = IMap(Energy2BPFormula());
	K.setParameter("mu",mu);

	g = IMap("var:u,s,alpha,I,v0,v1,eps;fun:u,s,alpha,v0,v1,0,eps;");

	A=A_;
	q=q_;
	h=h_;
	Kq=K(q)[0];

	to_R=IMap("var:alpha,J,r2,phi2,R2,PHI2,eps;fun:J*sin(alpha),J*cos(alpha),r2,phi2,R2,PHI2,eps;");
}

//   G(u,s,alpha,I,v0,v1,eps) =
//      (H(q+A(u,s,alpha,v0,v1,0,eps))-h,
//		 K(q+A(u,s,alpha,v0,v1,0,eps))-I-Kq)
//
//   G(y) =
//      (H(q+A*g(y))-h,
//		 K(q+A*g(y))-I-Kq).
IVector fromSectionCoordinates::G(const IVector &y) const
{
	IVector r(2);

	// y=(u,s,alpha,I,v0,v1,eps)
	interval I=y[3];
	
	r[0] = H(q+A*g(y))[0]-h;
	r[1] = K(q+A*g(y))[0]-I-Kq;
	return r;
}

IVector Y(const IVector &x,const IVector &v)
{
	IVector y(7);
	// v=(v0,v1)
	// x=(u,s,alpha,I,eps)
	// y=(u,s,alpha,I,v0,v1,eps)
	for(int i=0;i<4;i++) y[i]=x[i];
	y[4] = v[0];
	y[5] = v[1];
	y[6] = x[4];
	return y;
}
IVector fromSectionCoordinates::G(const IVector &x,const IVector &v) const
{
	return G(Y(x,v));
}

IMatrix fromSectionCoordinates::DG(const IVector &y) const
{
	IMatrix dG(2,7);
	
	IMatrix dH=H[q+A*g(y)]*A*g[y];
	IMatrix dK=K[q+A*g(y)]*A*g[y];
	
	// G(y) =
	//      (H(q+A*g(y))-h,
	//		 K(q+A*g(y))-I-Kq).
	// y=(u,s,alpha,I,v0,v1,eps)
	
	for(int i=0;i<7;i++) dG[0][i] = dH[0][i];
	for(int i=0;i<7;i++) dG[1][i] = dK[0][i];
	dG[1][3] = dK[0][3] - 1.0; // derivative with respect to I.
	return dG;
}

IMatrix fromSectionCoordinates::dGdx(const IVector &y) const
{
	IMatrix dG=DG(y);
	// x=(u,s,alpha,I,eps)
	// y=(u,s,alpha,I,v0,v1,eps)
	IMatrix D(2,5);
	for(int i=0;i<2;i++)
		for(int j=0;j<4;j++)
			D[i][j] = dG[i][j];
	D[0][4] = dG[0][6];
	D[1][4] = dG[1][6];
	return D;
}

IMatrix fromSectionCoordinates::dGdv(const IVector &y) const
{
	IMatrix dG=DG(y);
	// x=(u,s,alpha,I,eps)
	// y=(u,s,alpha,I,v0,v1,eps)
	IMatrix D(2,2);
	
	D[0][0] = dG[0][4]; D[0][1] = dG[0][5];
	D[1][0] = dG[1][4]; D[1][1] = dG[1][5];
	return D;
}

IMatrix fromSectionCoordinates::dGdv(const IVector &x,const IVector &v) const
{
	return dGdv(Y(x,v));
}

IMatrix fromSectionCoordinates::dGdx(const IVector &x,const IVector &v) const
{
	return dGdx(Y(x,v));
}

IVector fromSectionCoordinates::V(const IVector &x) const
{
	IVector v(2);
	IVector x0=midVector(x);
	// computing an initial approximation:
	for(int i=0;i<20;i++) v = midVector(v - gauss(dGdv(x0,v),G(x0,v)));
	
	// enlarging the initial approximation:
	for(int i=0;i<3;i++)
	{
		v = midVector(v) - gauss(dGdv(x,v),G(x,midVector(v)));
	}
	v = midVector(v) + 1.5*(v-midVector(v));

	// validation using Interval Newton method:
	IVector N = midVector(v) - gauss(dGdv(x,v),G(x,midVector(v)));
	if(subsetInterior(N,v)) return N;
	cout << "computation of the local change of coordinates failed. Aborting." << endl;
	abort();
	return IVector(2);
}

IMatrix fromSectionCoordinates::DV(const IVector &x) const
{
	IVector v=V(x);
	return -gaussInverseMatrix(dGdv(x,v))*dGdx(x,v);
}

IVector fromSectionCoordinates::omega(const IVector &x) const
{
	IVector v=V(x);
	IVector p(7);
	// x=(u,s,alpha,I,eps)
	// p=(u,s,alpha,v0,v1,0,eps)
	for(int i=0;i<3;i++) p[i] = x[i];
	p[3]=v[0];
	p[4]=v[1];
	p[5]=0.0;
	p[6]=x[4]; // eps
	return p;
}

IMatrix fromSectionCoordinates::Domega(const IVector &x) const
{
	IMatrix dv=DV(x);
	IMatrix D(7,5);
	// x=(u,s,alpha,I,eps)
	// p=(u,s,alpha,v0,v1,0,eps)
	for(int i=0;i<3;i++) D[i][i]=1.0;
	for(int i=0;i<2;i++)
		for(int j=0;j<5;j++)
			D[i+3][j] = dv[i][j];
	D[6][4]=1.0;
	return D;
}

///////////////////

fromLocalCoordinates::fromLocalCoordinates(IVector q,IMatrix A,interval h)
{
	psi=fromSectionCoordinates(q,A,h);
	B=IMatrix(5,5);
	for(int i=0;i<5;i++) B[i][i]=1.0;
}

fromLocalCoordinates::fromLocalCoordinates(IVector q,IMatrix A,IMatrix B_,interval h)
{
	psi=fromSectionCoordinates(q,A,h);
	B=B_;
}

///////////////

// toSectionCoordinates::toSectionCoordinates(IVector q_,IMatrix A)
// {
// 	q=q_;
// 	Ainv=gaussInverseMatrix(A);
// }

// IVector toSectionCoordinates::image(IVector p) const
// {
// 	p=Ainv*(p-q);

// 	IVector x(5);
// 	for(int i=0;i<4;i++) x[i]=p[i];
// 	x[4]=p[6]; // eps
// 	return x;
// }

// IMatrix toSectionCoordinates::derivative(const IVector &p) const
// {
// 	IMatrix DgInv(5,7);
// 	DgInv[0][0]=1;
// 	DgInv[1][1]=1;
// 	DgInv[2][2]=1; 
// 	DgInv[3][3]=1;
// 	DgInv[4][6]=1;
// 	return DgInv*Ainv;
// }

toSectionCoordinates::toSectionCoordinates(IVector q_,IMatrix A)
{
	q=q_;
	Ainv=gaussInverseMatrix(A);
	//K = IMap(sqrt2Energy2BPFormula());
	K = IMap(Energy2BPFormula());
	//K = IMap("par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:J;");
	K.setParameter("mu",mu);
	Kq=K(q)[0];
}

IVector toSectionCoordinates::image(IVector p) const
{
	IVector x(5);
	x[3]=K(p)[0]-Kq; // I
	p=Ainv*(p-q);
	for(int i=0;i<3;i++) x[i]=p[i]; // u,s,alpha
	x[4]=p[6]; // eps
	return x;
}

IMatrix toSectionCoordinates::derivative(const IVector &p) const
{
	IMatrix B(5,7);
	B[0][0]=1; // u
	B[1][1]=1; // s
	B[2][2]=1; // alpha
	B[4][6]=1; // eps
	IMatrix D=B*Ainv;
	IMatrix dK=K[p];
	for(int i=0;i<7;i++) D[3][i]=dK[0][i]; 
	return D;
}

////////////////////////////////

toLocalCoordinates::toLocalCoordinates(IVector q,IMatrix A)
{
	psiInv=toSectionCoordinates(q,A);
	Binv=IMatrix(5,5);
	for(int i=0;i<5;i++) Binv[i][i]=1.0;
}

toLocalCoordinates::toLocalCoordinates(IVector q,IMatrix A,IMatrix B)
{
	psiInv=toSectionCoordinates(q,A);
	Binv=gaussInverseMatrix(B);
}


