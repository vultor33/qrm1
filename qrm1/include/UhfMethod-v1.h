#ifndef UHFMETHOD_H
#define UHFMETHOD_H

#include <vector>
#include <fstream>
#include <string>

#include "Molecule.h"
#include "FourCenterMatrix.h"

// CUIDADO!!! ALGUNS PASSOS QUESTIONAVEIS

class UhfMethod
{
public:
	UhfMethod(){}
	~UhfMethod();
	void startUhfFockMatrix(int fock_matrix_size, Molecule & mol, FourCenterMatrix & four_center_, std::ofstream &escreve_);
	void buildFockMatrix(const std::vector< std::vector<double> > &coreFockMatrix);
	void build_first_density_matrix();
	void electronic_energy(std::vector< std::vector<double> > &coreFockMatrix);
	inline double getElectronicEnergy(){ return electronicEnergy; }

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
	void fock_matrix_diagonalization_and_density();
	void build_density_matrix(const std::vector< std::vector<double> > &eigenvectors, std::string typeMatrix);
	double electronicEnergy;

	//printing
	void printAll();
	void print_square_matrix(const std::vector< std::vector<double> > &matrix_to_be_printed, double unit_conversion);

	Molecule * pMol;
	FourCenterMatrix * pFourCenter_;
	std::ofstream * pEscreve_;

};

#endif