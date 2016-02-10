#ifndef OPTIMIZEWITHDLIB_H
#define OPTIMIZEWITHDLIB_H

#include "PrintAll.h"
#include "Molecule.h"
#include "ScfProcedure.h"

class OptimizeWithDlib
{
public:
	OptimizeWithDlib(){}
	~OptimizeWithDlib();

	void startOptimization();
	void optimizeH2(PrintAll &printLog_, Molecule &mol, bool printExcell);
	void optimizeH2p(PrintAll &printLog_, Molecule &mol, bool printExcell);
	void optimizeH3p(PrintAll &printLog_, Molecule &mol, bool printExcell);
	void optimizeH4(PrintAll &printLog_, Molecule &mol, bool printExcell);
	void optimizeH5p(PrintAll &printLog_, Molecule &mol, bool printExcell);

	void optimizen(
		PrintAll &printLog_, 
		Molecule &mol, 
		double finalEnergy,
		std::vector<int> &conections,
		std::vector<double> &startingPoint,
		ScfProcedure &scf_
		);
};

#endif