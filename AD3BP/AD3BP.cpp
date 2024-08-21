#include <iostream>
using namespace std;
#include <iomanip>
#include <omp.h>

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "constants.h"
#include "VectorField.h"
#include "readOrbits.h"
#include "localCoordinateChange.h"
#include "localMap.h"

interval part(interval x, int n, int i)
{
	return x.left()+i*(x.right()-x.left())/n+(x-x.left())/n;
}

interval mod2pi(interval t)
{
	if(t<2.0*interval::pi()) return t;
	return mod2pi(t-2.0*interval::pi());
}

interval strip()
{
	return interval(0.0,0.0825);
}

bool inStrip(interval alpha)
{
	return subsetInterior(mod2pi(alpha),strip());
}

IVector isCovering(const IVector &N,const IVector &fN,const IVector &fx0,const IMatrix &Df,bool &res)
{
	res=1;

	IVector M(5);
	M[0]=N[0];  // x
	M[1]=fN[1]; // y
	M[2]=fN[2]; // alpha
	M[3]=N[3]+Df[3][4]*N[4]; // I
	M[4]=N[4];  // eps

	IVector Nl=N;
	Nl[0]=Nl[0].left();
	IVector fNl = fx0 + Df*(Nl-midVector(N));
	if(not (fNl[0]<M[0]))
	{
		cout << "Problem with left image needed for a covering relation. " << endl;
		res=0;
	}

	IVector Nr=N;
	Nr[0]=Nr[0].right();
	IVector fNr = fx0 + Df*(Nr-midVector(N));
	if(not (fNr[0]>M[0]))
	{
		cout << "Problem with right image needed for a covering relation. " << endl;
		res=0;
	}

	return M;
}

IVector propagateCone(IVector cone,const IMatrix &Df,bool &res)
{
	res=1;
	cone=Df*cone;
	if(cone[0]<0.0) cone=-cone;
	if(not (cone[0]>0.0)) 
	{
		cout << "Problem with cone condition. " << endl;
		res=0;
	}
	for(int i=1;i<5;i++) cone[i]=cone[i]/cone[0].left();
	cone[0]=interval(1.0);
	return cone;
}

IVector initialHSet(interval alpha)
{
	IVector N(5);
	double r=pow(10.0,-9.0);
	N[0] = r*interval(-1,1); // x
	N[1] = r*interval(-1,1); // y
	N[2] = alpha;
	N[3] = pow(10.0,-11.0)*interval(0,1); // I
	N[4] = pow(10.0,-10.0)*interval(0,1); // eps
	return N;
}

interval initialHSetDiameter()
{
	return initialHSet(interval(0))[0].right();
}

IVector initialCone()
{
	IVector cone(5);
	cone[0]=interval(1.0);
	for(int i=1;i<5;i++) cone[i]=0.01*interval(-1,1);
	return cone;
}

interval initialConeSlope()
{
	return initialCone()[1].right();
}

bool coneCondition(IVector cone)
{
	IVector C=initialCone();
	for(int i=1;i<5;i++)
	{
		if(!(subsetInterior(cone[i],C[i]))) return 0;
	}
	return 1;
}

bool validateConnectingSequence(IMap &f,IOdeSolver &solver,interval alpha,interval &Ichange)
{
	vector<localMap> F=sequenceOfLocalMaps(f,solver);

	IVector N=initialHSet(alpha);
	IVector cone=initialCone();

	interval I0;
	interval dI;
	interval I(0.0);

	IVector fx0(5);
	IMatrix Df(5,5);

	int n=F.size();
	bool res=1;
	for(int i=0;i<n;i++)
	{	
		IVector fN=F[i](N,fx0,Df,dI);
		N = isCovering(N,fN,fx0,Df,res);
		if(res==0) return 0;

		cone = propagateCone(cone,Df,res);
		if(res==0) return 0;

		I=I+dI;
		// checking if we have returned to the strip after homoclinic excursion:
		if(i==97)
		{
			if(inStrip(N[2])) i=n; // exiting loop
		}
		// If we have not returned to the strip after 97th map, we make aditional turns
		// along the Lyapunov orbit and check again after the loop is finished.
		// By then we should be back in the strip.
	}

	if(!(inStrip(N[2]))) // checking if we returned to the strip
	{
		cout << "strip return failed." << endl;
		return 0; 
	}
	if(!(coneCondition(cone))) // checking if the propagated cone is tighter than the initial cone.
	{
		cout << "final cone condition failed." << endl;
		return 0;
	}
	if(!(I>0)) // We need to be increasing in I
	{
		cout << "Action condition failed." << endl;
		return 0;
	}
	if(!(subsetInterior(N[1],initialHSet(alpha)[1]))) // checking if we have covering (contraction along y-coordinate) for the final iterate of the map
	{
		cout << "final covering failed." << endl;
		return 0;
	}

	// When validateConnectingSequence() is initiaed for the first time,
	// then Ichange=0. In such case we set Ichange = I. Otherwise we take
	// Ichange = intervalHull(Ichange,I). In the way we set up the function
	// if we reach below condition we have I>0 so this way Ichange is computed only for
	// succesful runs and Ichange>0.
	if(Ichange>0)
	{
		Ichange = intervalHull(Ichange,I);
	}else
	{
		Ichange = I;
	}
	return 1;
}

void writeFinalResult(int flag,const vector<interval> &dI,double time)
{
	ofstream file("results/0_final_result.txt");
	file.precision(10);

	interval totalI=dI[0];
	int n=dI.size();
	for(int i=0;i<n;i++)
	{
		totalI=intervalHull(totalI,dI[i]);
	}
	
	if(flag==1)
	{
		file << "The proof was fully successful." << endl;
		file << "The total lower bound on the energy change is: " << totalI << endl;
	}else
	{
		file << "The proof was NOT fully succesful!" << endl;
		file << "Indexes of files which failed are in the sets failure_.txt. " << endl;
		file << "The total bound on the energy change for succesfull runs is: " << totalI << endl;
	}
	file << "Total computational time: " << time << endl;
	file << "Number of threads used: " << n << endl;
}

void timer(int count,int step,ofstream &file)
{
	if((count % step)==0) file << count << endl;
}

void singleRun()
{
	double itime, ftime, exec_time;
    itime = omp_get_wtime();

	IMap f(vectorFieldFormula());
	f.setParameter("mu",mu);
	IOdeSolver solver(f,TAYLOR_ORDER);
	interval dI;

	int M=200000;

	interval alpha=part(strip(),M,M/2)+2*initialConeSlope()*initialHSetDiameter()*interval(-1,1);
	
	if(validateConnectingSequence(f,solver,alpha,dI)==1)
	{
		cout << "the test for alpha="<< alpha << " was succesful." << endl;
	}

	ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("\n\nTime taken is %f", exec_time);
}

int main()
{
	clock_t start, end;
    double time;
    start = clock();
  	cout.precision(10);

	int M=200000;

	int N_of_threads=omp_get_max_threads();
	ofstream fileP("progress.txt");
	fileP << "Number of threads: " << N_of_threads << endl;

	vector<IMap*> f(N_of_threads); // vector field of the 3BP. We have these as a vector of objects, each object for a given processor. (To avoid potential clashes.)
	vector<IOdeSolver*> solver(N_of_threads); // these will be C^1 solvers, which use the vector fields f.
	vector<interval> dI(N_of_threads); // results for each thread
	vector<ofstream> fileS(N_of_threads);
	vector<ofstream> fileF(N_of_threads);
	ofstream fileE("results/errors.txt");

	for(int i=0;i<N_of_threads;i++)
	{
		f[i]=new IMap(vectorFieldFormula());
		f[i]->setParameter("mu",mu);
		solver[i] = new IOdeSolver(*f[i],TAYLOR_ORDER);
		fileS[i].open("results/succesful_"+to_string(i)+".txt");
		fileF[i].open("results/failure_"+to_string(i)+".txt");
	}

	int i, successFlag=1, count=0;

	#pragma omp parallel for private(i)
	for(i=0;i<M;i++)
	{
		interval alpha=part(strip(),M,i)+2*initialConeSlope()*initialHSetDiameter()*interval(-1,1);
			
		int id=omp_get_thread_num();
		try
		{
			if(validateConnectingSequence(*(f[id]),*(solver[id]),alpha,dI[id])==1)
			{
				fileS[id] << i << endl;
			}else
			{
				successFlag=0;
				fileF[id] << i << endl;
			}
		}catch(exception& e)
  		{
  			successFlag=0;
  			fileE << i << endl;
  			fileE << "Exception caught: "<< e.what() << endl << endl;
  		}
		timer(count,1000,fileP);
		count++;
	}
		
	for(int i=0;i<omp_get_max_threads();i++)
	{
		delete f[i];
		delete solver[i];
		fileS[i].close();
		fileF[i].close();
	}
	end = clock();
    time = (double(end) - double(start)) / CLOCKS_PER_SEC;

    writeFinalResult(successFlag,dI,time);
  	return 0;
} 
