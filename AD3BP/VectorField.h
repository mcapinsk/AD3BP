#ifndef VectorField_h
#define VectorField_h

// This function allows us to read the particular coefficients of the 
// vector field from files.
inline string formula(string fileName)
{
	string formula;
	ifstream file("1_VectorField/"+fileName+".txt");
	file >> formula;
	file.close();
	return formula;
}

// This returns the formula for the vector field of the full 3BP
// in the coordinates 
//    alpha,J,r2,phi2,R2,PHI2,eps
// REMARK: We add the parameter as last coordinate as a standing convention,
// throughout the entire code.
// 
// The formulae are derived in Mathematica in the
//    1_VectorField/MathematicaFormulaeComputation.nb
// file. This computes the vector field and stores the coefficients
// in .txt files. These formulae are then read by our program.
inline string vectorFieldFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("dalpha")+"," 
		+formula("dJ")+","
		+formula("dr2")+","
		+formula("dphi2")+","
		+formula("d_R2")+","
		+formula("d_PHI2")+",0.0;";
}

// This function returns the function which gives the value of the
// hamiltonian of the full 3BP computed for a point in the 
//     alpha,J,r2,phi2,R2,PHI2,eps
// coordinates.
inline string EnergyFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("H")+";";
}

// This function returns the function which gives the value of the
// Keplerian part K of the hamiltonian of the full 3BP computed,
//  for a point in the 
//     alpha,J,r2,phi2,R2,PHI2,eps
// coordinates.
inline string Energy2BPFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("Keps")+";";
}

#endif