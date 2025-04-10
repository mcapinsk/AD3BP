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
#include "readSections.h"
#include "localCoordinateChange.h"
#include "localMap.h"

// This finction subdivides an interval x into n equal parts,
// and returns its i-th part.
interval part(interval x, int n, int i)
{
	return x.left()+i*(x.right()-x.left())/n+(x-x.left())/n;
}

// This function computes a rigorous enclosure of an interval
// t module 2*PI.
interval mod2pi(interval t)
{
	if(t<2.0*interval::pi()) return t;
	return mod2pi(t-2.0*interval::pi());
}

// this function is used below for STEP 0. It returns the angle of our strip
interval strip()
{
	return interval(0.0,0.0825);
}
// This function is used for STEP 1. It gives the choice of the initial windows for 
// the connecting sequences, for the given choice of the angle interval.
IVector initialWindow(interval alpha)
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

// This chacks if the interval alpha lands inside or our strip:
bool inStrip(interval alpha)
{
	return subsetInterior(mod2pi(alpha),strip());
}

// This function validates a covering between two sets, N and M. 
// In the notations from the paper we validate that
//      f
//    N => M.
// The function in fact returns the set M, and ensures that it is 
// covered by N.
//
// We provide the function following:
//    N - te initial set for covering
//   fN - The interval enclosure of the image of the set N under the map f, i.e. f(N).
//  fx0 - The point x0 is the midpoint of N. The below function is given fx0 which is the 
//        enclosure for f(x0).
//   Df - is the interval encosure of the derivative of f on the set N, i.e. [Df(N)].
//  res - keeps track is the covering was successfuly valideted. If not,
//        this kept track of, and the computer assisted proof will report failure. 
//
// Comments: The function choosed M to be of the same size as N along the first coordinate.
// It also chooses the set M along the y,alpha so that we are sure that
// the image of N by f, i.e. f(N) when projected onto y,alpha coordinates
// lies in the y,alpha projection of M.
// We do the same for coordinate I, but here we use the maen value theorem.
// We know that for epsilon=0 the function f does not change I,
// so by the man value theorem 
//     \pi_I f(x,eps) = \pi_I (f(x,0.0)+df_deps(x,[0,eps])*eps)
//
// We then check that the images of the left and the right exit sets 
// of N are mapped out of M along the x-coordinate.
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

// Our cone is represented as 
//    cone = (1,a1,a2,a3,a4)
// where a1,a2,a3,a4 are closed intervals. This means that
//    C = {t*cone: t is real}
// The cone is propagated using the interval enclosure of 
// the derivative derivative. 
//
// The function returns another cone, so thet we are sure
// that the initial cone is returned into the returned one.
//
// The coefficients a1,a2,a3,a4 of the returned cone are computed 
// automatically.
//
// If the cone propagation fails, then "res" is set to zero and
// the computer assisted proof will report failure.
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

// This returns the size of the initial window N0 along the first coordinate.
// This is needed to choose an appropriate overlap of the windows in 
// the main part of the program.
interval initialWindowDiameter()
{
	return initialWindow(interval(0))[0].right();
}


// this is the cone we choose in N0. 
IVector initialCone()
{
	IVector cone(5);
	cone[0]=interval(1.0);
	for(int i=1;i<5;i++) cone[i]=0.01*interval(-1,1);
	return cone;
}

// This returns the initial cone slope.
// This is needed to choose an appropriate overlap of the windows in 
// the main part of the program.
interval initialConeSlope()
{
	return initialCone()[1].right();
}

// This function validates if a cone is sharper than the initialCone.
// This is needed to validate cone conditions at the final step of a
// connecting sequence.
bool coneCondition(IVector cone)
{
	IVector C=initialCone();
	for(int i=1;i<5;i++)
	{
		if(!(subsetInterior(cone[i],C[i]))) return 0;
	}
	return 1;
}

// This in many ways is the heart of the program. This function
// validates that we have the sequence of correct alignments
//   N0 => N1 => N2 => ... => Nk
// and that the local maps associated with this connecting sequence
// satisfy cone conditions.
// 
// What is also important is that the function measures the total change 
// in the coordinate I after the composition of the local maps involved 
// in the connecting sequence. What this means is that
//     I(f^k(x)) - I(x) \in eps*IchangeAlongSequence
// This value is then used to compute the GLOBAL Ichange, which is a bound
// that is valid for ALL connecting sequences.
// 
// Remark: For diffusion it is essential for us to ensure Ichange>0.
bool validateConnectingSequence(IOdeSolver &solver,interval alpha,interval &Ichange,bool vocal)
{
	vector<localMap> f=sequenceOfLocalMaps(solver);

	IVector N=initialWindow(alpha);
	IVector cone=initialCone();

	// here the variable dI measures the bound on the change of I for a single 
	// covering.
	interval dI;

	// This stores the total change of I along the connecting sequence
	interval IchangeAlongSequence(0.0);

	IVector fx0(5);
	IMatrix Df(5,5);

	int n=f.size();
	bool res=1;
	for(int i=0;i<n;i++)
	{	
		IVector fN=f[i](N,fx0,Df,dI);

		// STEP 2: Validation of the covering relation:
		N = isCovering(N,fN,fx0,Df,res);
		if(res==0) return 0;

		// STEP 2: Validation of the cone condition:
		cone = propagateCone(cone,Df,res);
		if(res==0) return 0;

		IchangeAlongSequence=IchangeAlongSequence+dI;

		// The user has the choice of running the whole proof or just
		// one connecting sequence. If we run just one sequence then below is 
		// displayed.
		if(vocal==1)
		{
			cout << "Validation of correct alignment of windows and cone conditions for the map f_" << i+1 << " OK. "<< endl;
			cout << "Bound on the change of I is " << IchangeAlongSequence << endl << endl;
		} 

		// We chack if we return to the strip. If we do, then there is no need to continue.
		if(i==97)
		{
			if(inStrip(N[2])) i=n; // exiting loop
		}
		// If we have not returned to the strip after 97th map, we make aditional turns
		// along the Lyapunov orbit and check again after the loop is finished.
		// By then we should be back in the strip.
	}
	// STEP 2: Validation of the final covering relation:
	if(!(subsetInterior(N[1],initialWindow(alpha)[1]))) // checking if we have covering (contraction along y-coordinate) for the final iterate of the map
		// this check ensures that we can in fact choose the final window N in the sequence to be of the same
		// size on the stable coordinate as the strip S.
	{
		cout << "final covering failed." << endl;
		return 0;
	}

	// STEP 3: checking if we returned to the strip
	if(!(inStrip(N[2]))) 
	{
		cout << "strip return failed." << endl;
		return 0; 
	}

	// STEP 4: checking if the propagated cone is tighter than the initial cone.
	if(!(coneCondition(cone))) 
	{
		cout << "final cone condition failed." << endl;
		return 0;
	}

	// STEP 5: The I is the bound on the energy change after passing 
	// through the connecting sequence. It is passed to Ichange, which computes
	// the global bound which holds for all the connecting sequences. 
	// We need to be increasing in I; this is checked below:
	if(!(IchangeAlongSequence>0)) 
	{
		cout << "Action condition failed." << endl;
		return 0;
	}
	// When validateConnectingSequence() is initiaed for the first time,
	// then Ichange=0. In such case we set Ichange = IchangeAlongSequence. Otherwise we take
	// Ichange = intervalHull(Ichange,IchangeAlongSequence). In the way we set up the function
	// if we reach below condition we have I>0 so this way Ichange is computed only for
	// succesful runs and Ichange>0.
	if(Ichange>0)
	{
		Ichange = intervalHull(Ichange,IchangeAlongSequence);
	}else
	{
		Ichange = IchangeAlongSequence;
	}
	return 1;
}

// The code is executed on multiple threads. On each thread we separately compute
// the Ichange along connecting sequences. Below function collects these bounds
// into a single "totalIchange".
void writeFinalResult(int flag,const vector<interval> &Ichange,double time)
{
	ofstream file("results/0_final_result.txt");
	file.precision(10);

	interval totalIchange=Ichange[0];
	int n=Ichange.size();
	for(int i=0;i<n;i++)
	{
		totalIchange=intervalHull(totalIchange,Ichange[i]);
	}
	
	if(flag==1)
	{
		file << "The proof was fully successful." << endl;
		file << "The total bound on the energy change is: " << totalIchange << endl;

		cout << "The proof was fully successful." << endl;
		cout << "The total bound on the energy change is: " << totalIchange << endl;
	}else
	{
		file << "The proof was NOT fully succesful!" << endl;
		file << "Indexes of files which failed are in the sets failure_.txt. " << endl;
		file << "The total bound on the energy change for succesfull runs is: " << totalIchange << endl;

		cout << "The proof was NOT fully succesful!" << endl;
		cout << "Indexes of files which failed are in the sets failure_.txt. " << endl;
		cout << "The total bound on the energy change for succesfull runs is: " << totalIchange << endl;
	}
	file << "Total computational time: " << time << endl;
	file << "Number of threads used: " << n << endl;

	cout << "Total computational time: " << time << endl;
	cout << "Number of threads used: " << n << endl;
}

// This function is executed when we call 
//    ./AD3BP k
// It validates a single connecting k-th sequence and displays the results.
void singleRun(int M,int k)
{
	IMap F(vectorFieldFormula());
	F.setParameter("mu",mu);
	IOdeSolver solver(F,TAYLOR_ORDER);
	interval dI;

	interval alpha=part(strip(),M,k)+2*initialConeSlope()*initialWindowDiameter()*interval(-1,1);
	
	cout << "Investigating " << k << "-th strip interval, where alpha=" << alpha << "." << endl << endl;
	if(validateConnectingSequence(solver,alpha,dI,1)==1)
	{
		cout << "The test for alpha="<< alpha << " was succesful." << endl;
		cout << "The final bound on the change of I is " << dI << endl << endl;
	}else
	{
		cout << "The required conditions failed! The proof is not succesful. "<< endl;
	}
}

int main(int argc, char* argv[])
{
	clock_t start, end;
	double time;
    start = clock();
  	cout.precision(10);

  	// Step 0. The choices of points and of the local coordinate changes are
  	// contained in folders: 
  	//   "00_points/midpoints.txt" this contains the points qi for i=0,...,122
  	//   "01_linear-coordinate-changes" The files "A.txt" and "B.txt" contain the matrices 
  	//   Ai and Bi, respectively. 

  	// This is the number of connecting sequences, which we will investigate.
  	int L=90000;

  	// This means that the program was called to validate a single selected 
  	// connecting sequence:
  	if(argc==2) 
  	{
  		int k = std::atoi(argv[1]);
  		if( (k>=0) and (k<L)) singleRun(L,k);
  		end = clock();
    	time = (double(end) - double(start)) / CLOCKS_PER_SEC;
  	}
  	// This means that the entire code will be executed:
  	else{
  		int N_of_threads=omp_get_max_threads();
		cout << "Number of threads: " << N_of_threads << endl;
		cout << "In total we need to validate " << L << " connecting sequences." << endl;
		cout << "We roport each consecutive 1000 connecting sequences which have been validated: " << endl;

		// For each thread we allocate a separate solver:
		vector<IMap*> F(N_of_threads); // vector field of the 3BP. We have these as a vector of objects, each object for a given processor. (To avoid potential clashes.)
		vector<IOdeSolver*> solver(N_of_threads); // these will be C^1 solvers, which use the vector fields f.
		vector<interval> Ichange(N_of_threads); // results for each thread
		// These files will store the results. 
		// fileS will store the succesful connecting sequences.
		// fileF will store the failed connecting sequences. (There should be none,
		// but we take into accound unforeseable events, and then if these files are not empty 
		// we will have information where things went wrong.)
		vector<ofstream> fileS(N_of_threads);
		vector<ofstream> fileF(N_of_threads);
		// We also have a separate file if erors are thrown along the computation.
		// This typically happens if in CAPD we divide by an interval containing zero or 
		// something is not computable. An error of this type also means the failure of the computer assisted proof.
		// (There should be no errors in a succesful run,
		// but we take into accound unforeseable events, and then if this file is not empty 
		// we will have information where things went wrong.)
		ofstream fileE("results/errors.txt");

		for(int i=0;i<N_of_threads;i++)
		{
			F[i]=new IMap(vectorFieldFormula());
			F[i]->setParameter("mu",mu);
			solver[i] = new IOdeSolver(*F[i],TAYLOR_ORDER);
			fileS[i].open("results/succesful_"+to_string(i)+".txt");
			fileF[i].open("results/failure_"+to_string(i)+".txt");
		}

		int l, successFlag=1, count=0;

		#pragma omp parallel for private(l)
		for(l=0;l<L;l++)
		{
			// STEP 1 Choice of the initial window.
			// The alpha is passed to validateConnectingSequence() where 
			// the window is creaded to be Bu x Bs x alpha x [0,10^{-11}]
			// This is done by calling N=initialWindow(alpha) inside of validateConnectingSequence(). 
			interval alpha=part(strip(),L,l)+2*initialConeSlope()*initialWindowDiameter()*interval(-1,1);
			
			int id=omp_get_thread_num();
			try
			{
				// STEP 2, STEP 3, STEP 4 and STEP 5 are performed in validateConnectingSequence().
				if(validateConnectingSequence(*(solver[id]),alpha,Ichange[id],0)==1)
				{
					fileS[id] << l << endl;
				}else
				{
					// The computer assisted proof failed! The program will continue running,
					// but the failure will be reported, and at the end the program will
					// write out that the computer assisted proof failed.
					successFlag=0;
					fileF[id] << l << endl;
				}
			}catch(exception& e)
  			{
  				// The computer assisted proof failed! The program will continue running,
				// but the failure will be reported, and at the end the program will
				// write out that the computer assisted proof failed.	
  				successFlag=0;
  				fileE << l << endl;
  				fileE << "Exception caught: "<< e.what() << endl << endl;
  			}
			count++;
			// We write out the progress:
			if((count % 1000)==0) cout << count << endl;
		}
		
		end = clock();
    	time = (double(end) - double(start)) / CLOCKS_PER_SEC;
		// STEP 6. The global bound on c and C from (66) is written into the result
		// data file.
		writeFinalResult(successFlag,Ichange,time);

		// The computer assisted proof is finished. All that remains is to clean up:
		for(int i=0;i<omp_get_max_threads();i++)
		{
			delete F[i];
			delete solver[i];
			fileS[i].close();
			fileF[i].close();
		}
  	}
  	cout << "computation time was: " << time << endl;
	
  	return 0;
} 
