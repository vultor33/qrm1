#ifndef SCFCYCLE_H
#define SCFCYCLE_H

#include <vector>
#include <fstream>
#include <string>

#include "MatrixDiagonalization.h"
#include "CoreMatrixCalculations.h"
#include "FourCenterMatrix.h"
#include "RhfMethod.h"
#include "UhfMethod.h"
#include "PrintAll.h"
#include "ScfMethod.h"

class ScfCycle
{
public:
	~ScfCycle();
	ScfCycle(Molecule &mol_in, PrintAll &logPrintLog_, std::string scfMethod_in = "RHF");
	bool SCF_step_i(int i_passo, int diagMethod = 0);
	double getFinalEnergy();
	double getIonizationPotential();
	bool qNotNan(){ return initializationSucess; }

private:
	int num_eletrons;
	int fock_matrix_size;
	bool converged;
	double core_core_repulsion();
	void scfCycleInitialization();
	bool initializationSucess;
	void printScfIteration(int iStep);
	std::string  scfMethod;
	std::vector< std::vector<double> > coreFockMatrix;
	double calculateRmsDensity();
	std::vector< std::vector<double> > oldDensity;
	double oldElectronicEnergy;
	bool checkDiis();
	double expq(double x, double q);

	Molecule &mol; //reference of outside Molecule object
	FourCenterMatrix four_center_;
	MatrixDiagonalization matrDiag_;//(remove) just to symmetrize core matrix
	RhfMethod rhfMatrix_;
	UhfMethod uhfMatrix_;
	ScfMethod scfMeth_;
	PrintAll * pPrintLog_;


};

#endif

