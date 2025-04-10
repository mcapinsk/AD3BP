#include "readSections.h"

///////////////////////////////////////////////
// FUNCTION: q0()
//
// This is our choice of the point q_0. It plays a special role,
// since our strip is positioned on the section at this point.
// We therefore write out this particular vector explicitly in the code.
IVector q0()
{
	IVector q(7);
	q[0]=interval(0); // alpha
	q[1]=interval(1)/interval(10); // J
	q[2]=0.951431;	// r2
	q[3]=interval::pi(); // phi2
	q[4]=interval(0); // R2
	q[5]=0.9708302420604081; // PHI2
	q[6]=interval(0); // eps
	return q;
}

///////////////////////////////////////////////
// FUNCTION: readOrbit()
//
// As described in the .h file.
vector<IVector> readOrbit()
{
	ifstream inFile("0_sections/q.txt");
	int n=123;
	vector<IVector> q(n);
	
	// This is an auxiliary vector used to read from the file:
	IVector v(7); 
	// This is an auxiliary vector used to read from the file:
	double x;	

	for(int k=0;k<n;k++)
	{
		v[0]=0.0; // alpha
		v[1]=interval(1)/interval(10); // J
		inFile >> x; v[2]=x; // r2
		inFile >> x; v[3]=x; // phi2
		inFile >> x; v[4]=x; // R2
		inFile >> x; v[5]=x; // PHI2
		v[6]=0.0; // eps
		q[k]=v;
	}

	q[0]=q0();
	q[98]=q0();
	q[n-1]=q0();
	
	return q;
}

///////////////////////////////////////////////
// FUNCTION: scalarProduct()
//
// REMARK: The role of this function is self-explanatory, but one issue is
// worthwhile commenting. The function returns an interval. This interval
// is ensured to contain the true scalar product between all the vectors 
// contained in the IVectors x,y.  
interval scalarProduct(IVector x,IVector y)
{
	interval s(0);
	int n=x.dimension();
	for(int i=0;i<n;i++) s=s+x[i]*y[i];
	return s;
}

///////////////////////////////////////////////
// FUNCTION: orthogonalPart()
//
// This function takes two vectors as arguments, v and n, and returns
// the projection of v onto the orthogonal space to n.
// 
// REMARK: The function is suitable for computer assisted proofs in the sense
// that if the function returns an IVector, say X, then we are sure that this
// IVector contains a true vector x, for which 
//   scalarProduct(x,n) = 0.
IVector orthogonalPart(IVector v,IVector n)
{
	return v - (scalarProduct(v,n)/scalarProduct(n,n))*n;
}

///////////////////////////////////////////////
// FUNCTION: column()
// 
// This function returns k-th colimn of a matrix X.
IVector column(const IMatrix &A,int k)
{
	int n=A.numberOfRows();
	IVector v(n);
	for(int i=0;i<n;i++) v[i]=A[i][k];
	return v;
}

/////////////////////////////////////////////////
// FUNCTION: makeOrthogonalToLastColumn()
//
// This function modifies a given matrix A, so that 
// its k-th column becomes orthogonal to the column with index 5.
//
// REMARK: The reason why we choose index 5 will be explained in the 
// comments for the next routine.
void makeOrthogonalToFifthColumn(IMatrix &A,int k)
{
	IVector v=orthogonalPart(column(A,k),column(A,5));
	int n=v.dimension();
	for(int i=0;i<n;i++) A[i][k]=v[i];
}

/////////////////////////////////////////////////
// FUNCTION: makeOrthogonalToLastColumn()
//
// This function modifies a given matrix A, so that 
// columns its columns indexed by 0,1,2,3,4 become orthogonal
// to the column with index 5. As outlined in the paper, the columns
// of the matrix A consist of vectors v0,...,v5, where 
//   v0 - is a vector aligned towards the unstable direction
//   v1 - is a vector aligned towards the stable direction
//   v2 - is (alpha,J,r2,phi2,R2,PHI2,epsilon) = (1,0,...,0)
//   v3 - is (alpha,J,r2,phi2,R2,PHI2,epsilon) = (0,1,0,...,0)
//   v4 - is (alpha,J,r2,phi2,R2,PHI2,epsilon) = (0,0,GradientH0,0)
//   v5 - is (alpha,J,r2,phi2,R2,PHI2,epsilon) = (0,0,F0,0)
//   v6 - is the additional collumn corresponding to the
//        fact that we add parameter as the last variable.
//        So v6 is (alpha,J,r2,phi2,R2,PHI2,epsilon) = (0,...,0,1).
// All this is not aparent in the code. In fact, the code is
// indifferent to our choice of the matrix A in the sense that
// if we choose some different matrix and all computations go through,
// then the computer assisted proof will remain valid.
// What is essential is that the columns v0,...v4,v6 are orthogonal to v5.
// This is essential because the way we compute section-to-section maps
// (which are in CAPD named as Poincare maps) relies on the fact that the 
// section is orthogonal to the chosen vector. In our applications we will
// choose the sections to be orthogonal to v5.
void makeOrthogonalToFifthColumn(IMatrix &A)
{
	for(int i=0;i<5;i++) makeOrthogonalToFifthColumn(A,i);
	makeOrthogonalToFifthColumn(A,6);	
}

///////////////////////////////////////////////
// FUNCTION: A0()
//
// This is our choice of the matrix A_0. It plays a special role,
// since our strip is positioned on the section Sigma0 which is defined by this matrix.
// We therefore write out this particular matrix explicitly in the code.
IMatrix A0()
{
	IMatrix A(7,7);

	A[0][0]=0.98586255112737;	A[0][1]=0.9858597924844;	A[0][2]=interval(1);	A[0][3]=interval(0);	A[0][4]=interval(0);		A[0][5]=interval(0);		A[0][6]=interval(0);
	A[1][0]=interval(0);		A[1][1]=interval(0);		A[1][2]=interval(0);	A[1][3]=interval(1);	A[1][4]=interval(0);		A[1][5]=interval(0);		A[1][6]=interval(0);
	A[2][0]=-0.075305747051793;	A[2][1]=0.07530553633099;	A[2][2]=interval(0);	A[2][3]=interval(0);	A[2][4]=-0.079682214206858;	A[2][5]=interval(0);		A[2][6]=interval(0);
	A[3][0]=0.09224658009712;	A[3][1]=0.092246321972443;	A[3][2]=interval(0);	A[3][3]=interval(0);	A[3][4]=interval(0);		A[3][5]=0.072478762739209;	A[3][6]=interval(0);
	A[4][0]=-0.08390728168075;	A[4][1]=-0.083907046890831;	A[4][2]=interval(0);	A[4][3]=interval(0);	A[4][4]=interval(0);		A[4][5]=0.079682214206858;	A[4][6]=interval(0);
	A[5][0]=-0.082790164191965; A[5][1]=0.082789932528175;	A[5][2]=interval(0);	A[5][3]=interval(0);	A[5][4]=0.072478762739209;	A[5][5]=interval(0);		A[5][6]=interval(0);
	A[6][0]=interval(0);		A[6][1]=interval(0);		A[6][2]=interval(0);	A[6][3]=interval(0);	A[6][4]=interval(0);		A[6][5]=interval(0);		A[6][6]=interval(1);
	
	makeOrthogonalToFifthColumn(A);
	return A;
}

/////////////////////////////////////////////////
// FUNCTION: getLinearChanges()
// 
// This routine returns the sequence of matrices A_i, which define
// the local surfaces of sections. For each i=0,...,122 a section
// will be chosen to be orthogonal to the column of A_i whose index is 5.
// We read the matrices from the file 
//    01_linear-coordinate-changes/01_linear-coordinate-changes/A.txt
// Such matrices were computed numerically, so we ensure that all the columns
// with indeces 0,1,2,3,4,6 are orthogonal to the column with index 5.
//
// We also make sure that A_0 = A_98 = A_122 = A0().
vector<IMatrix> getLinearChanges()
{
	ifstream file("0_sections/A.txt");
	
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
		makeOrthogonalToFifthColumn(D);
		A[k]=D;
	}
	A[0]=A0();
	A[98]=A0();
	A[n-1]=A0();
	return A;
}

///////////////////////////////////////////////
// FUNCTION: B0()
//
// This is our choice of the matrix B_0. It plays a special role,
// since our strip is positioned on the section Sigma0. 
// The matrix B0 introduced local coordinates on this section.
// This is why we write it out explicitly in the code.
//
// REMARK: We can choose any matrix of the required form, which means that we can
// freely choose the coefficients
//    B[0][3], B[0][4],
//    B[1][3], B[1][4].
// As long as the computer program validates the needed conditions, this choice is arbitrary.
// This particular choice aligns well the dynamics expressed in the local coordinates.
// Our particular choice of these coefficients comes from a careful, non-rigorous
// numerical investigation of the system.
// 
// The above comments are true also for the matrices B[k] which are read in the
// below routine getLocalLinearChanges().
//
// The last, 5-th column of the matrix B0 (and of the matrices B[k] below)
// corresponds to the coordinate epsilon, which is added as the last variable.
// What this means is that the last column of the matrix stores the vectoe
// w_0 from (34) from the paper.
//
// REMARK:
// The theory states that we also have freedom of choice of the coefficients
//   B[0][2],
//   B[1][2].
// We set these to zero.
IMatrix B0()
{
	IMatrix B(5,5);
	B[0][0]=interval(1);	B[0][1]=interval(0);	B[0][2]=interval(0); 	B[0][3]=35.320864504461;	B[0][4]=-0.2482024853181;
	B[1][0]=interval(0);	B[1][1]=interval(1);	B[1][2]=interval(0);	B[1][3]=-35.320913922206;	B[1][4]=-0.26775285797094;
	B[2][0]=interval(0);	B[2][1]=interval(0);	B[2][2]=interval(1);	B[2][3]=interval(0);		B[2][4]=interval(0);
	B[3][0]=interval(0);	B[3][1]=interval(0);	B[3][2]=interval(0);	B[3][3]=interval(1);		B[3][4]=interval(0);
	B[4][0]=interval(0);	B[4][1]=interval(0);	B[4][2]=interval(0);	B[4][3]=interval(0);		B[4][4]=interval(1);
	return B;
}

/////////////////////////////////////////////////
// FUNCTION: getLocalLinearChanges()
//
// This function returns the matrices which are used to define the local
// coordinates on the sections. We make sure that the matrices are of required form.
// Namely, they have to be equal to one on the diagonal, and the only freedom of choice we have
// is for the coefficients
//    B[k][0][3], B[k][0][4],
//    B[k][1][3], B[k][1][4].
//
// The last, 5-th column of the matrices B[k] corresponds to the coordinate epsilon, 
// which is added as the last variable.
// What this means is that the last column of the matrix B[k] stores the vectoe
// w_k from (34) from the paper.
//
// We make sure that B_0 = B_98 = B_122 = B0().
//
// REMARK:
// The theory states that we also have freedom of choice of the coefficients
//   B[k][0][2],
//   B[k][1][2].
// We set these to zero.
vector<IMatrix> getLocalLinearChanges()
{
	ifstream file("0_sections/B.txt");
	
	int n=123;
	vector<IMatrix> B(n);
	
	IMatrix D(5,5);
	double x;
	for(int k=0;k<n;k++)
	{
		// we read a matrix from the file
		for(int i=0;i<5;i++)
		{
			for(int j=0;j<5;j++)
			{
				file >> x;
				D[i][j]=x;
			}
		}

		B[k]=IMatrix(5,5); // this creates a 5x5 matrix filled with zeros.
		// We make sure that B[k] has ones on the diagonal.
		for(int i=0;i<5;i++) B[k][i][i]=interval(1);

		// We choose these coefficients from the data file:
		B[k][0][3] = D[0][3]; B[k][0][4] = D[0][4];
		B[k][1][3] = D[1][3]; B[k][1][4] = D[1][4];  
	}
	B[0]=B0();
	B[98]=B0();
	B[n-1]=B0();
	return B;
}