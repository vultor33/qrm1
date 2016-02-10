#ifndef SCFMETHOD_H
#define SCFMETHOD_H

#include <vector>
#include <fstream>
#include <string>

#include "RhfMethod.h"
#include "UhfMethod.h"
#include "Molecule.h"
#include "FourCenterMatrix.h"


class ScfMethod
{
public:
	ScfMethod(){ sucessfulStart = false; }
	~ScfMethod();

	void startScfMethod(std::string method_in, int fock_matrix_size, Molecule & mol, FourCenterMatrix & four_center_, PrintAll * pPrintLog_in_);
	void scfMethodCycle(int diagMethod, std::vector< std::vector<double> > &coreFockMatrix);
	void setDiis(bool newDiis);
	void safeDeletePointers();

	double getElectronicEnergy();
	std::vector< std::vector<double> > getDensity();
	bool getDiis();
	double getIonizationPotential();

	//safe delete criar um delete so com tudo

private:
	bool sucessfulStart;
	std::string method;
	RhfMethod rhfMatrix_;
	UhfMethod uhfMatrix_;
	
};

#endif


