#ifndef readSections_h
#define readSections_h

#include "capd/capdlib.h"
#include <iostream>
using namespace capd;
using namespace std;

// This routine returns the sequence of points q_0,...,q_122 at which we position
// the surfaces of section. They are precomputed and stored in 
//    00_points/midpoints.txt
//
// REMARK: The routine makes sure that we choose
//   q_0=q0()
//   q_98=q0()
//   q_122=q0()
// where the vector q0() is given in the first routine of the readOrbits.cpp file.
// 
// REMARK: Our code is arranged so that the perturbation paremeter epsilon is 
// included as the last variable. So, the vectors q_i are seven dimensional. We 
// set the last, epsilon-coordinate to be zero.
//
// REMARK: We ensure that all the points have J=1/10 and alpha=0. 
// 
// REMARK: The points are in the coordinates (alpha,J,r2,phi2,R2,PHI2,epsilon).
vector<IVector> readOrbit();

// This routine returns the sequence of matrices A_0,...,A_122 which define 
// surfaces of section. 
vector<IMatrix> getLinearChanges();

// This routine returns the sequence of matrices B_0,...,B_122 which define
// the local coordinates on the surfaces of section. 
vector<IMatrix> getLocalLinearChanges();

// This routine returns the point q0 at which the section Sigma0 is positioned.
// (The section Sigma0 plays a special role in our construction, since this is wher
// the strip is positioned.)
// We make this accessible to other parts of the program, since it is used
// in localCoordinateChange.cpp to compute the energy kappa0=K(q0). 
// It is essential that all local coordinate sections use the same kappa0.
// Passing the same q0 is our way of ensuring this.
IVector q0();

#endif