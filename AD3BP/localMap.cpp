#include "localMap.h"
#include "readOrbits.h"
#include "VectorField.h"
#include "constants.h"

localMap::localMap(IVector q0,IVector q1,IMatrix A0,IMatrix A1,IMatrix B0,IMatrix B1,interval h_,IMap &f_,IOdeSolver &solver_)
{
	q.resize(2); A.resize(2); B.resize(2);
	q[0]=q0; q[1]=q1;
	A[0]=A0; A[1]=A1;
	B[0]=B0; B[1]=B1;
	
	h=h_;
	
	psi   =fromLocalCoordinates(q[0],A[0],B[0],h);
	psiInv=toLocalCoordinates(q[1],A[1],B[1]);
	
	IVector w(7);
	for(int i=0;i<7;i++) w[i]=A[1][i][5];
	
	section = new IAffineSection(q[1],w);
	f = &f_;
	solver = &solver_;//new IOdeSolver(*f,TAYLOR_ORDER);
	P = new IPoincareMap(*solver,*section);
}

localMap::~localMap()
{
	delete P;
	delete section;
}

localMap& localMap::operator=(const localMap& other)
{
	if (this != &other)
	{
		q=other.get_q();
		A=other.get_A();
		B=other.get_B();
		h=other.get_h();

		f = other.f; 
		solver = other.solver; 
		
		psi   =fromLocalCoordinates(q[0],A[0],B[0],h);
		psiInv=toLocalCoordinates(q[1],A[1],B[1]);
	
		IVector w(7);
		for(int i=0;i<7;i++) w[i]=A[1][i][5];
		section = new IAffineSection(q[1],w);

		P = new IPoincareMap(*solver,*section);
		K = new IMap(Energy2BPFormula());
		K->setParameter("mu",mu);
	}
	return *this;
}

IVector localMap::image(const IVector &x,IVector &fx0,IMatrix &Df,interval &dI) const
{
	IVector x0=midVector(x);

	C0Rect2Set Q0=psi.C0Set(x0);
	C1Rect2Set Q1=psi.C1Set(x);

	IMatrix DP(7,7);
	IVector y=(*P)(Q1,DP);

	DP = P->computeDP(y,DP);
	Df = psiInv[y]*DP*psi[x];
	dI = ((*K)[y]*DP-(*K)[psi(x)])[0][6];
	// in fact, we could compute the interval enclosure of the change 
	// of energy as Df[3][4], but this would be far less accurate 
	// than the above computation of dI. (The above computed dI will be
	// a subset of Df[3][4].)
	Df[3][4] = dI; 

	try
	{
		fx0 = psiInv((*P)(Q0));
	}catch(exception& e)
  	{
  		C1Rect2Set C1Q0=psi.C1Set(x0);
  		fx0 = psiInv((*P)(C1Q0));
  	}

	return fx0 + Df*(x-x0);
}

//////////////////////

vector<localMap> sequenceOfLocalMaps(IMap &f,IOdeSolver &solver)
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
	
	interval h=H(q[0])[0];
	vector<localMap> F(n-1);
	
	for(int i=0;i<n-1;i++) F[i]=localMap(q[i],q[i+1],A[i],A[i+1],B[i],B[i+1],h,f,solver);
	return F;
}
