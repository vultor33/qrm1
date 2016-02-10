#ifndef RHFMETHOD_H
#define RHFMETHOD_H

#include <vector>
#include <fstream>

#include "Molecule.h"
#include "FourCenterMatrix.h"
#include "MatrixDiagonalization.h"
#include "PrintAll.h"
#include "DiisProcedure.h"

class RhfMethod
{
public:
	RhfMethod(){}
	~RhfMethod();
	void startRhfFockMatrix(int fock_matrix_size, Molecule & mol, FourCenterMatrix & four_center_, PrintAll * pPrintLog_in_);
	void buildFockMatrix(const std::vector< std::vector<double> > &coreFockMatrix);
	void electronic_energy(std::vector< std::vector<double> > &coreFockMatrix);
	void build_first_density_matrix();
	inline double getElectronicEnergy(){ return electronicEnergy; }
	inline void setDiagonalizationMethod(int diagMethod_in){ diagonalizationMethod = diagMethod_in; }
	std::vector< std::vector<double> > getDensity();
	double getIonizationPotential();
	void deletePointers();

//	void safeDiisPointerDelete();
	inline bool getDiis(){ return pDiis_->diisOn; }
	inline void setDiis(bool newDiis){ pDiis_->diisOn = newDiis; }

private:
	std::vector< std::vector<double> > fockMatrix;
	std::vector< std::vector<double> > densityMatrix;

	// Build fock matrix
	int fockMatrixSize;
	std::vector < int > AtomFockMap;
	std::vector < int > BaseFockMap;
	double calculateFockMatrix(int iFock, int jFock);
	double couloumbExchangeMiMi(int iFock);
	double couloumbExchangeMiLambda(int iFock,int jFock);
	double couloumbExchangeMiNi(int iFock,int jFock);

	//Diagonalization and density
	int num_electrons;
	int diagonalizationMethod;
	void fock_matrix_diagonalization_and_density();
	void build_density_matrix(const std::vector< std::vector<double> > &eigenvectors);
	double electronicEnergy;

	//objects
	Molecule * pMol;
	FourCenterMatrix * pFourCenter_;
	MatrixDiagonalization matrDiag_;
	PrintAll * pPrintLog_;
	DiisProcedure * pDiis_;

};

#endif
