#ifndef FITNESSOPT_H
#define FITNESSOPT_H

#include <vector>
#include <fstream>
#include <string>

#include "ScfProcedure.h"
#include "Molecule.h"

class FitnessOpt
{
public:
	FitnessOpt(){}
	~FitnessOpt();

	void printFitness(ScfProcedure &scf_, bool printExell, int model, std::string errorName);

private:
	double refEnergy;
	double refIonizPot;
	std::vector<double> refDistances;
	void setSystems(int model);
	double calculateDistanceError(Molecule *pMol);
	void printExcellDistances(std::ofstream &excell_, Molecule *pMOl);

};

#endif