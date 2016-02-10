#ifndef UHFMETHOD_H
#define UHFMETHOD_H

#include <vector>
#include <fstream>
#include <string>

#include "Molecule.h"
#include "FourCenterMatrix.h"
#include "MatrixDiagonalization.h"
#include "PrintAll.h"
#include "DiisProcedure.h"

// CUIDADO!!! ALGUNS PASSOS QUESTIONAVEIS

class UhfMethod
{
public:
	UhfMethod(){}
	~UhfMethod();
	void startUhfFockMatrix(int fock_matrix_size, Molecule & mol, FourCenterMatrix & four_center_, PrintAll * pPrintLog_in_);
	void buildFockMatrix(const std::vector< std::vector<double> > &coreFockMatrix);
	void build_first_density_matrix();
	void electronic_energy(std::vector< std::vector<double> > &coreFockMatrix);
	inline double getElectronicEnergy(){ return electronicEnergy; }
	inline void setDiagonalizationMethod(int diagMethod_in){ diagonalizationMethod = diagMethod_in; }
	void deletePointers();

	void setDiis(bool newDiis);
	bool getDiis();
	std::vector< std::vector<double> > getDensity();
	double getIonizationPotential();

private:
	std::vector< std::vector<double> > alphaFockMatrix;
	std::vector< std::vector<double> > betaFockMatrix;
	std::vector< std::vector<double> > alphaDensityMatrix;
	std::vector< std::vector<double> > betaDensityMatrix;

	// Build fock matrix
	int fockMatrixSize;
	std::vector < int > AtomFockMap;
	std::vector < int > BaseFockMap;
	double calculateFockMatrix(int iFock, int jFock, std::string typeMatrix);
	double couloumbExchangeMiMi(int iFock, std::string typeMatrix);
	double couloumbExchangeMiLambda(int iFock, int jFock, std::string typeMatrix);
	double couloumbExchangeMiNi(int iFock, int jFock, std::string typeMatrix);

	//Diagonalization and density
	int num_electrons;
	int alphaElectrons;
	int betaElectrons;
	int diagonalizationMethod;
	void fock_matrix_diagonalization_and_density();
	void build_density_matrix(const std::vector< std::vector<double> > &eigenvectors, std::string typeMatrix);
	double electronicEnergy;


	Molecule * pMol;
	FourCenterMatrix * pFourCenter_;
	MatrixDiagonalization matrDiagAlpha_;
	MatrixDiagonalization matrDiagBeta_;
	PrintAll * pPrintLog_;
	DiisProcedure * pDiisAlpha_;
	DiisProcedure * pDiisBeta_;

};

#endif