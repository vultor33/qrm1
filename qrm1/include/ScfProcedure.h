#ifndef SCFPROCEDURE_H
#define SCFPROCEDURE_H

#include <string>

#include "ScfCycle.h"
#include "Molecule.h"

class ScfProcedure
{
public:
	~ScfProcedure();
	ScfProcedure(){}

	void startScfProcedure(Molecule * mol_in, PrintAll * printLog_in_);

	void doScfProcedure(std::string inputName);
	double doScfProcedure(std::vector<double> &xyzCoord);
	double doScfProcedure(std::vector<CoordXYZ> &xyzCoord);
	//double

	inline Molecule *getPmol(){ return pMol; }
	inline double getIonizationPotential(){ return ionizationPotential; }
	inline double getFinalEnergy(){ return finalEnergy; }

private:
	Molecule * pMol;
	PrintAll * pPrintLog_;

	double ionizationPotential;
	double heatOfFormation;
	double finalEnergy;
};

#endif

