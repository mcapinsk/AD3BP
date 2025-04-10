#include "localMap.h"
#include "readSections.h"
#include "VectorField.h"
#include "constants.h"

localMap::localMap(IVector q0,IVector q1,IMatrix A0,IMatrix A1,IMatrix B0,IMatrix B1,IOdeSolver &solver_)
{
	q.resize(2); A.resize(2); B.resize(2);
	q[0]=q0; q[1]=q1;
	A[0]=A0; A[1]=A1;
	B[0]=B0; B[1]=B1;
	
	psi   =fromLocalCoordinates(q[0],A[0],B[0]);
	psiInv=toLocalCoordinates(q[1],A[1],B[1]);
	
	// We take the section1 ro be at q[1] and orthogonal to the
	// fifth column of A[1]. The vector w is chosen as the fifth column.
	IVector w(7);
	for(int i=0;i<7;i++) w[i]=A[1][i][5];
	
	section1 = new IAffineSection(q[1],w);
	solver = &solver_; 
	P = new IPoincareMap(*solver,*section1);
}

localMap::~localMap()
{
	delete P;
	delete section1;
}

localMap& localMap::operator=(const localMap& other)
{
	if (this != &other)
	{
		q=other.get_q();
		A=other.get_A();
		B=other.get_B();

		solver = other.solver; 
		
		psi   =fromLocalCoordinates(q[0],A[0],B[0]);
		psiInv=toLocalCoordinates(q[1],A[1],B[1]);
	
		IVector w(7);
		for(int i=0;i<7;i++) w[i]=A[1][i][5];
		section1 = new IAffineSection(q[1],w);

		P = new IPoincareMap(*solver,*section1);
		K = new IMap(Energy2BPFormula());
		K->setParameter("mu",mu);
	}
	return *this;
}

IVector localMap::image(const IVector &X0,IVector &fx0,IMatrix &Df,interval &dI) const
{
	IVector x0=midVector(X0);

	// CAPD performs computations more accurately if it integrates on
	// parallelogram representations C0Rect2Set and C1Rect2Set. The C0Rect2Set
	// is sed for C^0 computations and C1Rect2Set for C^1 computations.
	// The functions psi.C0Set(x0) and psi.C01et(x0) return parallelogram bounds
	// on the image of the box X0 passed to alpha,J,r2,phi2,R2,PHI2,eps coordinates.
	// See localCoordinateChange.h/cpp for more comments.
	C0Rect2Set Q0=psi.C0Set(x0);
	C1Rect2Set Q1=psi.C1Set(X0);

	// The DP will be the derivative of the map P to section1.
	IMatrix DP(7,7);
	// this computes the bound y = P(psi(X0))
	// as well as the monodromy matrix DP along the flow.
	IVector y=(*P)(Q1,DP); 
	// The monodromy matrix DP is recomputed to the derivative of P.
	DP = P->computeDP(y,DP);
	// Below we compute Df(X0). Recall that the operators [] compute derivatives
	// in the classes from localCoordinateChange.h/cpp.
	Df = psiInv[y]*DP*psi[X0];
	// Here we use (64) from the paper to compute:
	dI = ((*K)[y]*DP-(*K)[psi(X0)])[0][6];

	// It is sufficient for us to compute C^0 bound on f(x0).
	// We have observed that in some cases the integrator has problems
	// in performing the computation, hence we try doing this,
	try
	{
		fx0 = psiInv((*P)(Q0));
	}catch(exception& e)
	// and if computation of f(x0) is not validated using the C^0 solver,
	// we use C^1 computation, which is slower, but more accurate.
  	{
  		C1Rect2Set C1Q0=psi.C1Set(x0);
  		fx0 = psiInv((*P)(C1Q0));
  	}

  	// We use the mean value theorem to return f(X0):
	return fx0 + Df*(X0-x0);
}

//////////////////////

vector<localMap> sequenceOfLocalMaps(IOdeSolver &solver)
{
	IMap H(EnergyFormula());
	H.setParameter("mu",mu);
	
	vector<IVector> q=readOrbit();
	int n=q.size();

	vector<IMatrix> A=getLinearChanges();
	vector<IMatrix> B=getLocalLinearChanges();
	A[n-1]=A[0];
	B[n-1]=B[0];
	q[n-1]=q[0];
	
	vector<localMap> f(n-1);
	
	for(int i=0;i<n-1;i++) f[i]=localMap(q[i],q[i+1],A[i],A[i+1],B[i],B[i+1],solver);
	return f;
}
