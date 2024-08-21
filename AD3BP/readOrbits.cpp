#include "readOrbits.h"

IVector q0()
{
	IVector q(7);
	q[0]=interval(0);
	q[1]=interval(1)/interval(10);
	q[2]=0.951431;	
	q[3]=interval::pi();
	q[4]=interval(0);
	q[5]=0.9708302420604081;
	q[6]=interval(0);
	return q;
}

vector<IVector> readOrbit()
{
	ifstream inFile("00_points/midpoints.txt");
	int n=123;
	
	double x;
	vector<IVector> V(n);
	
	IVector v(7);
	for(int k=0;k<n;k++)
	{
		v[0]=0.0; // alpha
		v[1]=interval(1)/interval(10); // J
		inFile >> x; v[2]=x; // r2
		inFile >> x; v[3]=x; // phi2
		inFile >> x; v[4]=x; // R2
		inFile >> x; v[5]=x; // PHI2
		v[6]=0.0; // eps
		V[k]=v;
	}

	V[0]=q0();
	V[98]=q0();
	V[n-1]=q0();
	
	return V;
}

///////////////////////////////////////////

interval scalarProduct(IVector x,IVector y)
{
	interval s(0);
	int n=x.dimension();
	for(int i=0;i<n;i++) s=s+x[i]*y[i];
	return s;
}

IVector orthogonalPart(IVector v,IVector n)
{
	return v - (scalarProduct(v,n)/scalarProduct(n,n))*n;
}

IVector column(const IMatrix &A,int k)
{
	int n=A.numberOfRows();
	IVector v(n);
	for(int i=0;i<n;i++) v[i]=A[i][k];
	return v;
}

void makeOrthogonalToLastColumn(IMatrix &A,int k)
{
	IVector v=orthogonalPart(column(A,k),column(A,5));
	int n=v.dimension();
	for(int i=0;i<n;i++) A[i][k]=v[i];
}

void makeOrthogonalToLastColumn(IMatrix &A)
{
	for(int i=0;i<5;i++) makeOrthogonalToLastColumn(A,i);
}

vector<IMatrix> getLinearChanges()
{
	ifstream file("01_linear-coordinate-changes/A.txt");
	
	int n=123;
	vector<IMatrix> A(n);
	
	IMatrix D(7,7);
	double x;
	for(int k=0;k<n;k++)
	{
		for(int i=0;i<7;i++)
		{
			for(int j=0;j<7;j++)
			{
				file >> x;
				D[i][j]=x;
			}
		}
		makeOrthogonalToLastColumn(D);
		A[k]=D;
	}
	A[98]=A[0];
	A[n-1]=A[0];
	return A;
}

/////////////////////////////////////////

vector<IMatrix> getLocalLinearChanges()
{
	ifstream file("01_linear-coordinate-changes/B.txt");
	
	int n=123;
	vector<IMatrix> B(n);
	
	IMatrix D(5,5);
	double x;
	for(int k=0;k<n;k++)
	{
		for(int i=0;i<5;i++)
		{
			for(int j=0;j<5;j++)
			{
				file >> x;
				D[i][j]=x;
			}
		}
		B[k]=D;
	}
	B[98]=B[0];
	B[n-1]=B[0];
	return B;
}