#ifndef VectorField_h
#define VectorField_h

inline string formula(string fileName)
{
	string formula;
	ifstream file("02_VectorField/"+fileName+".txt");
	file >> formula;
	file.close();
	return formula;
}

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

inline string inverseVectorFieldFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		"-("+formula("dalpha")+"),"
		"-("+formula("dJ")+"),"
		"-("+formula("dr2")+"),"
		"-("+formula("dphi2")+"),"
		"-("+formula("d_R2")+"),"
		"-("+formula("d_PHI2")+"),0.0;";
}

inline string EnergyFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("H")+";";
}

inline string EnergyR3BPFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("Heps")+";";
}

inline string sqrt2Energy2BPFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		"(2.0*("+formula("Keps")+"))^0.5;";
}

inline string Energy2BPFormula()
{
	return "par:mu;var:alpha,J,r2,phi2,R2,PHI2,eps;fun:"
		+formula("Keps")+";";
}

#endif