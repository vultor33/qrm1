#include "RhfMethod.h"

#include <fstream>
#include <iomanip>

#include "Molecule.h"
#include "Params.h"
#include "FourCenterMatrix.h"
#include "MatrixDiagonalization.h"
#include "PrintAll.h"

using namespace std;

RhfMethod::~RhfMethod(){}

void RhfMethod::startRhfFockMatrix(int fock_matrix_size, 
	Molecule & mol, 
	FourCenterMatrix & four_center_, 
	PrintAll * pPrintLog_in_)
{
	pMol = &mol;
	pFourCenter_ = &four_center_;
	pPrintLog_ = pPrintLog_in_;
	num_electrons = mol.number_of_electrons;
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

	pDiis_ = new DiisProcedure(pPrintLog_, fockMatrix, densityMatrix);
}

void RhfMethod::buildFockMatrix(const vector< vector<double> > &coreFockMatrix)
{
	fockMatrix = coreFockMatrix;

	for (int iFock = 0; iFock < fockMatrixSize; iFock++)
	{
		for (int jFock = iFock; jFock < fockMatrixSize; jFock++)
		{
			fockMatrix[iFock][jFock] += calculateFockMatrix(iFock, jFock);
		}
	}
	
	pDiis_->diisCalculation();

	fock_matrix_diagonalization_and_density();

	pPrintLog_->printScfMatrix("fock", fockMatrix, Params::hartree_ev);
	pPrintLog_->printScfMatrix("density", densityMatrix);
}

double RhfMethod::calculateFockMatrix(int iFock, int jFock)
{
	if (iFock == jFock)
	{
		return(couloumbExchangeMiMi(iFock));
	}
	else
	{
		if (AtomFockMap[iFock] == AtomFockMap[jFock])
		{
			return(couloumbExchangeMiNi(iFock, jFock));
		}
		else
		{
			return(couloumbExchangeMiLambda(iFock, jFock));
		}
	}
}

double RhfMethod::couloumbExchangeMiMi(int iFock)
{
	double auxsoma = 0;
	int i_densidade;
	int j_densidade;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];

	for (int ni = 0; ni < Params::get_int(pMol->atom_name[A], "base_number"); ni++)
	{
		auxsoma += densityMatrix[pMol->i_fock_base_line[A][ni]][pMol->i_fock_base_line[A][ni]] *
			(pFourCenter_->get_four_center(A, A, mi, mi, ni, ni) -
			0.5e0*pFourCenter_->get_four_center(A, A, mi, ni, mi, ni));
	}

	for (int B = 0; B < pMol->number_of_atoms; B++)
	{
		if (B != A)
		{
			for (int lambda = 0; lambda < Params::get_int(pMol->atom_name[B], "base_number"); lambda++)
			{
				for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
				{
					i_densidade = pMol->i_fock_base_line[B][lambda];
					j_densidade = pMol->i_fock_base_line[B][sigma];
					auxsoma += densityMatrix[i_densidade][j_densidade] *
						pFourCenter_->get_four_center(A, B, mi, mi, lambda, sigma);
				}
			}
		}
	}
	return auxsoma;
}

double RhfMethod::couloumbExchangeMiLambda(int iFock, int jFock)
{
	//(int A, int B, int mi, int lambda)
	// ni e o lambda nesse caso
	double auxsoma = 0.0e0;
	int i_densidade;
	int j_densidade;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];
	int B = AtomFockMap[jFock];
	int lambda = BaseFockMap[jFock];

	for (int ni = 0; ni < Params::get_int(pMol->atom_name[A], "base_number"); ni++)
	{
		for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
		{
			i_densidade = pMol->i_fock_base_line[A][ni];
			j_densidade = pMol->i_fock_base_line[B][sigma];
			auxsoma += densityMatrix[i_densidade][j_densidade] *
				pFourCenter_->get_four_center(A, B, mi, ni, lambda, sigma, true);
		}
	}

	auxsoma *= -0.5e0;
	return auxsoma;
}

double RhfMethod::couloumbExchangeMiNi(int iFock, int jFock)
{

	double auxsoma = 0.0e0;
	int A = AtomFockMap[iFock];
	int mi = BaseFockMap[iFock];
	int ni = BaseFockMap[jFock];
	int i_densidade = pMol->i_fock_base_line[A][mi];
	int j_densidade = pMol->i_fock_base_line[A][ni];

	auxsoma += 0.5e0*densityMatrix[i_densidade][j_densidade] *
		(3.0e0 * pFourCenter_->get_four_center(A, A, mi, ni, mi, ni) -
		pFourCenter_->get_four_center(A, A, mi, mi, ni, ni));

	for (int B = 0; B < pMol->number_of_atoms; B++)
	{
		if (B != A)
		{
			for (int lambda = 0; lambda < Params::get_int(pMol->atom_name[B], "base_number"); lambda++)
			{
				for (int sigma = 0; sigma < Params::get_int(pMol->atom_name[B], "base_number"); sigma++)
				{
					int i_densidade = pMol->i_fock_base_line[B][lambda];
					int j_densidade = pMol->i_fock_base_line[B][sigma];

					auxsoma += densityMatrix[i_densidade][j_densidade] *
						pFourCenter_->get_four_center(A, B, mi, ni, lambda, sigma);
				}
			}
		}
	}

	return auxsoma;
}

void RhfMethod::build_density_matrix(const vector< vector<double> > &eigenvectors)
{
	/*
	fredmudar
	vector< vector<double> > testEigen;
	testEigen = eigenvectors;
	testEigen[0][0] = -0.93921125;
	testEigen[0][1] = -0.27521184;
	testEigen[0][2] = -3.10346112E-02;
	testEigen[0][3] = -0.20292248;
	testEigen[1][0] = -0.23194227;
	testEigen[1][1] = 0.15568359;
	testEigen[1][2] = 0.56567973;
	testEigen[1][3] = 0.77586848;
	testEigen[2][0] = -0.13018401;
	testEigen[2][1] = 0.75904691;
	testEigen[2][2] = 0.40905359;
	testEigen[2][3] = -0.48946396;
	testEigen[3][0] = -0.21711092;
	testEigen[3][1] = 0.56909472;
	testEigen[3][2] = -0.71534497;
	testEigen[3][3] = 0.34245530;
	*/

	for (int i = 0; i < fockMatrixSize; i++)
	{
		for (int j = 0; j < fockMatrixSize; j++)
		{
			densityMatrix[i][j] = 0.0e0;
			for (int alfa = 0; alfa < (num_electrons / 2); alfa++)
			{
				densityMatrix[i][j] += eigenvectors[i][alfa] * eigenvectors[j][alfa];
				//densityMatrix[i][j] += testEigen[i][alfa] * testEigen[j][alfa];
			}
			densityMatrix[i][j] *= 2;
		}
	}
}

void RhfMethod::fock_matrix_diagonalization_and_density()
{
	matrDiag_.diagonalization(fockMatrix, diagonalizationMethod, num_electrons / 2);
	build_density_matrix(matrDiag_.getEigenvectors());
}

void RhfMethod::electronic_energy(vector< vector<double> > &coreFockMatrix)
{
	electronicEnergy = 0.0e0;
	int fockSize = fockMatrix.size();
	//szabo pag 150
	for (int mi = 0; mi < fockSize ; mi++)
	{
		for (int ni = 0; ni < fockSize; ni++)
		{
			electronicEnergy += densityMatrix[mi][ni] * (fockMatrix[mi][ni] + coreFockMatrix[mi][ni]);
		}
	}
	electronicEnergy *= 0.5e0;
}

void RhfMethod::build_first_density_matrix()
{
	densityMatrix.resize(fockMatrixSize);
	for (int def_densidade = 0; def_densidade < fockMatrixSize; def_densidade++)
	{
		densityMatrix[def_densidade].resize(fockMatrixSize);
	}

	int charge = pMol->getCharge();
	double chargeWeight = 1.0e0;
	if(charge != 0 )
		chargeWeight = (double)pMol->number_of_electrons / (double)(pMol->number_of_electrons + charge);

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
					densityMatrix[i][j_base] = (double)0.0e0;
				}
				else
				{
					densityMatrix[i][i] = chargeWeight * ((double)Params::get_int(pMol->atom_name[B], "number_of_electrons")) / ((double)Params::get_int(pMol->atom_name[B], "base_number"));
				}
			}
		}
	}
	pPrintLog_->printScfMatrix("first density", densityMatrix);
}

vector< vector<double> > RhfMethod::getDensity()
{
	return densityMatrix;
}

double RhfMethod::getIonizationPotential()
{
	int ionizationPos = -1 + num_electrons / 2;
	return -matrDiag_.getEigenvalueI(ionizationPos);
}

void RhfMethod::deletePointers()
{
	delete pDiis_;
}


/*
void RhfMethod::safeDiisPointerDelete()
{
vector< vector<double> > dummyFock;
vector< vector<double> > dummyDensity;
pDiis_ = new DiisProcedure(NULL, dummyFock, dummyDensity);
}
*/




