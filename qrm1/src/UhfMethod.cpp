#include "UhfMethod.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>

#include "Molecule.h"
#include "Params.h"
#include "FourCenterMatrix.h"
#include "MatrixDiagonalization.h"

using namespace std;

UhfMethod::~UhfMethod(){}

void UhfMethod::startUhfFockMatrix(int fock_matrix_size, 
	Molecule & mol, 
	FourCenterMatrix & four_center_, 
	PrintAll * pPrintLog_in_)
{
	pMol = &mol;
	pFourCenter_ = &four_center_;
	pPrintLog_ = pPrintLog_in_;
	num_electrons = mol.number_of_electrons;
	int charge = mol.getCharge();
	num_electrons -= charge;
	if ((num_electrons % 2) == 1)
	{
		alphaElectrons = 1 + (num_electrons - 1) / 2;
		betaElectrons = (num_electrons - 1) / 2;
	}
	else
	{
		alphaElectrons = num_electrons / 2;
		betaElectrons = num_electrons / 2;
	}


	diagonalizationMethod = 0;

	fockMatrixSize = fock_matrix_size;
	AtomFockMap.resize(fock_matrix_size);
	BaseFockMap.resize(fock_matrix_size);

	int iFock = 0;
	for (int A = 0; A < pMol->number_of_atoms; A++)
	{
		int n_bases_A = Params::get_int(pMol->atom_name[A], "base_number");
		for (int mi = 0; mi < n_bases_A; mi++)
		{
			AtomFockMap[iFock] = A;
			BaseFockMap[iFock] = mi;
			iFock++;
		}
	}

	pDiisAlpha_ = new DiisProcedure(pPrintLog_, alphaFockMatrix, alphaDensityMatrix);
	pDiisBeta_ = new DiisProcedure(pPrintLog_, betaFockMatrix, betaDensityMatrix);
}

void UhfMethod::buildFockMatrix(const vector< vector<double> > &coreFockMatrix)
{
	alphaFockMatrix = coreFockMatrix;
	betaFockMatrix = coreFockMatrix;

	for (int iFock = 0; iFock < fockMatrixSize; iFock++)
	{
		for (int jFock = iFock; jFock < fockMatrixSize; jFock++)
		{
			alphaFockMatrix[iFock][jFock] += calculateFockMatrix(iFock, jFock,"alpha");
			betaFockMatrix[iFock][jFock] += calculateFockMatrix(iFock, jFock, "beta");
		}
	}

	pDiisAlpha_->diisCalculation();
	pDiisBeta_->diisCalculation();

	fock_matrix_diagonalization_and_density();

	//impressao fredmudar
	pPrintLog_->printScfMatrix("fock", alphaFockMatrix, Params::hartree_ev);
	pPrintLog_->printScfMatrix("density", alphaDensityMatrix);
	pPrintLog_->printScfMatrix("fock", betaFockMatrix, Params::hartree_ev);
	pPrintLog_->printScfMatrix("density", betaDensityMatrix);

}

double UhfMethod::calculateFockMatrix(int iFock, int jFock, string typeMatrix)
{
	if (iFock == jFock)
	{
		return(couloumbExchangeMiMi(iFock,typeMatrix));
	}
	else
	{
		if (AtomFockMap[iFock] == AtomFockMap[jFock])
		{
			return(couloumbExchangeMiNi(iFock, jFock, typeMatrix));
		}
		else
		{
			return(couloumbExchangeMiLambda(iFock, jFock, typeMatrix));
		}
	}
}

double UhfMethod::couloumbExchangeMiMi(int iFock, string typeMatrix)
{
	double auxsoma = 0;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];

	for (int ni = 0; ni < Params::get_int(pMol->atom_name[A], "base_number"); ni++)
	{
		if (typeMatrix == "alpha")
		{
			auxsoma += (
				(alphaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][mi]] +
				betaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][mi]])*
				pFourCenter_->get_four_center(A, A, mi, mi, ni, ni)
				) -
				alphaDensityMatrix[pMol->i_fock_base_line[A][ni]][pMol->i_fock_base_line[A][ni]] *
				pFourCenter_->get_four_center(A, A, mi, ni, mi, ni);
		}
		else
		{
			auxsoma += (
				(alphaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][mi]] +
				betaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][mi]])*
				pFourCenter_->get_four_center(A, A, mi, mi, ni, ni)
				) -
				betaDensityMatrix[pMol->i_fock_base_line[A][ni]][pMol->i_fock_base_line[A][ni]] *
				pFourCenter_->get_four_center(A, A, mi, ni, mi, ni);
		}

	}

	for (int B = 0; B < pMol->number_of_atoms; B++)
	{
		if (B != A)
		{
			for (int lambda = 0; lambda < Params::get_int(pMol->atom_name[B], "base_number"); lambda++)
			{
				for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
				{
					auxsoma +=
						(alphaDensityMatrix[pMol->i_fock_base_line[B][lambda]][pMol->i_fock_base_line[B][sigma]] +
						betaDensityMatrix[pMol->i_fock_base_line[B][lambda]][pMol->i_fock_base_line[B][sigma]])*
						pFourCenter_->get_four_center(A, B, mi, mi, lambda, sigma);
				}
			}
		}
	}
	return auxsoma;
}

double UhfMethod::couloumbExchangeMiLambda(int iFock, int jFock, string typeMatrix)
{
	//(int A, int B, int mi, int lambda)
	// ni e o lambda nesse caso
	double auxsoma = 0.0e0;
	int i_densidade;
	int j_densidade;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];
	int B = AtomFockMap[jFock];
	int ni = BaseFockMap[jFock];

	for (int lambda = 0; lambda < Params::get_int(pMol->atom_name[A], "base_number"); lambda++)
	{
		for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
		{
			i_densidade = pMol->i_fock_base_line[A][lambda];
			j_densidade = pMol->i_fock_base_line[B][sigma];
			if (typeMatrix == "alpha")
			{
				auxsoma += alphaDensityMatrix[i_densidade][j_densidade] *
					pFourCenter_->get_four_center(A, B, mi, ni, lambda, sigma);
			}
			else
			{
				auxsoma += betaDensityMatrix[i_densidade][j_densidade] *
					pFourCenter_->get_four_center(A, B, mi, ni, lambda, sigma);
			}
		}
	}
	auxsoma *= -1.0e0;
	return auxsoma;
}

double UhfMethod::couloumbExchangeMiNi(int iFock, int jFock, string typeMatrix)
{

	double auxsoma = 0.0e0;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];
	int ni = BaseFockMap[jFock];
	int i_densidade = pMol->i_fock_base_line[A][mi];
	int j_densidade = pMol->i_fock_base_line[A][ni];



	if (typeMatrix == "alpha")
	{
		auxsoma += (
			2.0e0*(alphaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]] +
			betaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]])*
			pFourCenter_->get_four_center(A, A, mi, ni, mi, ni)
			) -
			alphaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]] *
			(
			pFourCenter_->get_four_center(A, A, mi, ni, mi, ni) +
			pFourCenter_->get_four_center(A, A, mi, mi, ni, ni)
			);
	}
	else
	{
		auxsoma += (
			2.0e0*(alphaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]] +
			betaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]])*
			pFourCenter_->get_four_center(A, A, mi, ni, mi, ni)
			) -
			betaDensityMatrix[pMol->i_fock_base_line[A][mi]][pMol->i_fock_base_line[A][ni]] *
			(
			pFourCenter_->get_four_center(A, A, mi, ni, mi, ni) +
			pFourCenter_->get_four_center(A, A, mi, mi, ni, ni)
			);
	}

	for (int B = 0; B < pMol->number_of_atoms; B++)
	{
		if (B != A)
		{
			for (int lambda = 0; lambda < Params::get_int(pMol->atom_name[B], "base_number"); lambda++)
			{
				for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
				{
					auxsoma +=
						(alphaDensityMatrix[pMol->i_fock_base_line[B][lambda]][pMol->i_fock_base_line[B][sigma]] +
						betaDensityMatrix[pMol->i_fock_base_line[B][lambda]][pMol->i_fock_base_line[B][sigma]])*
						pFourCenter_->get_four_center(A, B, mi, ni, lambda, sigma);
				}
			}
		}
	}
	return auxsoma;
}

void UhfMethod::build_density_matrix(const vector< vector<double> > &eigenvectors, string typeMatrix)
{
	for (int i = 0; i < fockMatrixSize; i++)
	{
		for (int j = 0; j < fockMatrixSize; j++)
		{
			if (typeMatrix == "alpha")
			{
				alphaDensityMatrix[i][j] = 0.0e0;
				for (int alfa = 0; alfa < alphaElectrons; alfa++)
				{
					alphaDensityMatrix[i][j] += eigenvectors[i][alfa] * eigenvectors[j][alfa];
				}
			}
			else
			{
				betaDensityMatrix[i][j] = 0.0e0;
				for (int beta = 0; beta < betaElectrons; beta++)
				{
					betaDensityMatrix[i][j] += eigenvectors[i][beta] * eigenvectors[j][beta];
				}
			}

		}
	}
}

void UhfMethod::fock_matrix_diagonalization_and_density()
{
	matrDiagAlpha_.diagonalization(alphaFockMatrix, diagonalizationMethod, alphaElectrons);
	build_density_matrix(matrDiagAlpha_.getEigenvectors(), "alpha");
	matrDiagBeta_.diagonalization(betaFockMatrix, diagonalizationMethod, betaElectrons);
	build_density_matrix(matrDiagBeta_.getEigenvectors(), "beta");
}

void UhfMethod::electronic_energy(vector< vector<double> > &coreFockMatrix)
{
	electronicEnergy = 0.0e0;
	for (int mi = 0; mi < fockMatrixSize; mi++)
	{
		for (int ni = 0; ni < fockMatrixSize; ni++)
		{
			electronicEnergy += 
				((alphaDensityMatrix[mi][ni] + betaDensityMatrix[mi][ni])*
				coreFockMatrix[mi][ni])+
				alphaDensityMatrix[mi][ni] * alphaFockMatrix[mi][ni]+
				betaDensityMatrix[mi][ni] * betaFockMatrix[mi][ni];
		}
	}
	electronicEnergy *= 0.5e0;
}

void UhfMethod::build_first_density_matrix()
{
	alphaDensityMatrix.resize(fockMatrixSize);
	betaDensityMatrix.resize(fockMatrixSize);
	for (int setDensity = 0; setDensity < fockMatrixSize; setDensity++)
	{
		alphaDensityMatrix[setDensity].resize(fockMatrixSize);
		betaDensityMatrix[setDensity].resize(fockMatrixSize);
	}
	int charge = pMol->getCharge();

//	double chargeWeight = 1.0e0;
//	if (charge != 0)
//		chargeWeight = pMol->number_of_electrons / (pMol->number_of_electrons - charge);
	
	double alphaFraction = (double)alphaElectrons / ((double)num_electrons);
	double betaFraction = (double)betaElectrons / ((double)num_electrons);

	int j_base = 0;
	for (int i = 0; i < fockMatrixSize; i++)
	{
		for (int B = 0; B < pMol->number_of_atoms; B++)
		{
			for (int j = 0; j < Params::get_int(pMol->atom_name[B], "base_number"); j++)
			{
				j_base = pMol->i_fock_base_line[B][j];
				if (i != j_base)
				{
					alphaDensityMatrix[i][j_base] = (double)0.0e0;
					betaDensityMatrix[i][j_base] = (double)0.0e0;
				}
				else
				{
					alphaDensityMatrix[i][i] = alphaFraction*((double)Params::get_int(pMol->atom_name[B], "number_of_electrons")) / ((double)Params::get_int(pMol->atom_name[B], "base_number"));
					betaDensityMatrix[i][i] =  betaFraction*((double)Params::get_int(pMol->atom_name[B], "number_of_electrons")) / ((double)Params::get_int(pMol->atom_name[B], "base_number"));
				}
			}
		}
	}

	pPrintLog_->printScfMatrix("first density", alphaDensityMatrix);
	pPrintLog_->printScfMatrix("first density", betaDensityMatrix);
}

vector< vector<double> > UhfMethod::getDensity()
{
	int size = alphaDensityMatrix.size();
	vector< vector<double> > density;
	density.resize(size);
	for (int i = 0; i < size; i++)
		density[i].resize(size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			density[i][j] = alphaDensityMatrix[i][j] + betaDensityMatrix[i][j];
		}
	}
	return density;
}

double UhfMethod::getIonizationPotential()
{
	int ionizationPos = -1 + alphaElectrons;
	return -matrDiagAlpha_.getEigenvalueI(ionizationPos);
}

bool UhfMethod::getDiis()
{
	return pDiisAlpha_->diisOn;
}

void UhfMethod::setDiis(bool newDiis)
{
	pDiisAlpha_->diisOn = newDiis;
	pDiisBeta_->diisOn = newDiis;
}

void UhfMethod::deletePointers()
{
	delete pDiisAlpha_;
	delete pDiisBeta_;
}





