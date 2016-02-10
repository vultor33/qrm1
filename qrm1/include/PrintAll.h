#ifndef PRINTALL_H
#define PRINTALL_H

#include <vector>
#include <fstream>
#include <string>

#include "Coordstructs.h"

class PrintAll
{
public:
	PrintAll(std::string outputName, int debugLevel = 0);
	~PrintAll();

	void printString(std::string anyString);

	void printOpening();

	void printScfMatrix(
		std::string matrixType,
		const std::vector< std::vector<double> > &entryMatrix,
		double unitConversion = 1.0e0);
	void printScfMatrix(
		std::string matrixType,
		const std::vector<double> &entryMatrix,
		double unitConversion = 1.0e0);


	void printMolecule(
		const std::vector<CoordZMAT> & zmatMolecule,
		const std::vector<CoordXYZ> & xyzMolecule,
		int nAtoms
		);

	void printNonZeroIntegrals(
		int nAtoms,
		std::vector< std::vector<diatomicFourCenter> > &matrixIntegrals,
		std::vector< std::string > & atomName
		);

	void printIteration(int iStep, double energy, double energyVariation = 0.0e0, double densityRms = 0.0e0);
	void printScfHeader(int flag = 0);
	void printFinalEnergy(double elecEnergy, double coreEnergy);
	void printEndOfScf(bool converged);
	void printStartDiis();
	void printErrorDiis(double maxError);
	void printEndDiis();

private:
	std::string outputName;
	std::ofstream logOutput_;
	//debug
	void debugOptions(int debugLevel);
	bool conditionFirst;
	bool conditionMatrixIterations;
	bool conditionXyzCoordinates;
	bool conditionDiis;
	bool checkPrintConditions(std::string type);

	void printMatrix(const std::vector< std::vector<double> > &entryMatrix, double unitConversion);
	void printOnXyzFile(const std::vector<CoordXYZ> & xyzMolecule, int nAtoms);
	void printXyzOnLogOutput(const std::vector<CoordXYZ> & xyzMolecule, int nAtoms);
	void printZmatOnLogOutput(const std::vector<CoordZMAT> & zmatMolecule, int nAtoms);


};


#endif