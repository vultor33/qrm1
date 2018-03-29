#include "ScfCycle.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <string>

#include "MatrixDiagonalization.h"
#include "CoreMatrixCalculations.h"
#include "Params.h"
#include "FourCenterMatrix.h"
#include "RhfMethod.h"
#include "PrintAll.h"

using namespace std;

ScfCycle::ScfCycle(
	Molecule &mol_in, 
	PrintAll &printLog_,
	string scfMethod_in)
:mol(mol_in)
{
	pPrintLog_ = &printLog_;
	scfMethod = scfMethod_in;
	converged = false;
	scfCycleInitialization();
}

ScfCycle::~ScfCycle(){}

bool ScfCycle::SCF_step_i(int iStep, int diagMethod)
{
	scfMeth_.scfMethodCycle(diagMethod, coreFockMatrix);
	printScfIteration(iStep);
	return converged;
}

double ScfCycle::getFinalEnergy()
{
	double coreEnergy = core_core_repulsion();
	pPrintLog_->printFinalEnergy(Params::hartree_ev*oldElectronicEnergy, coreEnergy, Params::hartree_ev*getIonizationPotential());
	return Params::hartree_ev*oldElectronicEnergy + coreEnergy;
}

double ScfCycle::core_core_repulsion()
{
	double core_core_repulsion = 0.0e0;
	double RAB;

	double alfA;
	double a1A, b1A, c1A, a2A, b2A, c2A, a3A, b3A, c3A, a4A, b4A, c4A;
	double F_A;
	
	double alfB;
	double a1B, b1B, c1B, a2B, b2B, c2B, a3B, b3B, c3B, a4B, b4B, c4B;
	double F_B;

	double za_zb;
	double gamass;

	for (int A = 0; A < (mol.number_of_atoms-1); A++)
	{
		for (int B = (A+1); B < mol.number_of_atoms; B++)
		{
			//Rab de ua para angstrons
			RAB =  Params::bohr_angs * mol.distancia_i_j[A][B];
			alfA = Params::get_double(mol.atom_name[A],"alfacore");
			a1A = Params::get_double(mol.atom_name[A], "a1core");
			b1A = Params::get_double(mol.atom_name[A], "b1core");
			c1A = Params::get_double(mol.atom_name[A], "c1core");
			a2A = Params::get_double(mol.atom_name[A], "a2core");
			b2A = Params::get_double(mol.atom_name[A], "b2core");
			c2A = Params::get_double(mol.atom_name[A], "c2core");
			a3A = Params::get_double(mol.atom_name[A], "a3core");
			b3A = Params::get_double(mol.atom_name[A], "b3core");
			c3A = Params::get_double(mol.atom_name[A], "c3core");
			a4A = Params::get_double(mol.atom_name[A], "a4core");
			b4A = Params::get_double(mol.atom_name[A], "b4core");
			c4A = Params::get_double(mol.atom_name[A], "c4core");

			/*Agora a parte de B
			AM1: A New General Purpose Quantum Mechanical
			Molecular Model’ - Michael Dewar, 1984
			*/

			alfB = Params::get_double(mol.atom_name[B], "alfacore");
			a1B = Params::get_double(mol.atom_name[B], "a1core");
			b1B = Params::get_double(mol.atom_name[B], "b1core");
			c1B = Params::get_double(mol.atom_name[B], "c1core");
			a2B = Params::get_double(mol.atom_name[B], "a2core");
			b2B = Params::get_double(mol.atom_name[B], "b2core");
			c2B = Params::get_double(mol.atom_name[B], "c2core");
			a3B = Params::get_double(mol.atom_name[B], "a3core");
			b3B = Params::get_double(mol.atom_name[B], "b3core");
			c3B = Params::get_double(mol.atom_name[B], "c3core");
			a4B = Params::get_double(mol.atom_name[B], "a4core");
			b4B = Params::get_double(mol.atom_name[B], "b4core");
			c4B = Params::get_double(mol.atom_name[B], "c4core");


			//sohidrogeniomudar
			double qgauss = Params::get_double(mol.atom_name[A], "qgauss");
			string coreMethod = Params::method;

			if(coreMethod == "skewRM1")
			{
				double d1A = Params::get_double(mol.atom_name[A], "skew");
				double d1B = Params::get_double(mol.atom_name[B], "skew");
				F_A = gaussSkew(RAB,a1A,b1A,c1A,d1A);
				F_B = gaussSkew(RAB,a1B,b1B,c1B,d1B);
			}
			else
			{
				if (
					(coreMethod[0] == 'q')
					)
				{
					F_A = a1A*expq((-b1A*(RAB - c1A)*(RAB - c1A)), qgauss);
					F_B = a1B*expq((-b1B*(RAB - c1B)*(RAB - c1B)), qgauss);
				}
				else
				{
					F_A = a1A*exp(-b1A*(RAB - c1A)*(RAB - c1A));
					F_B = a1B*exp(-b1B*(RAB - c1B)*(RAB - c1B));
				}
					if (
					((coreMethod != "RM1-1g") && (coreMethod != "qRM1-1g"))
					&& (coreMethod[0] != 'q')
					)
				{
					F_A += a2A*exp(-b2A*(RAB - c2A)*(RAB - c2A));
					F_B += a2B*exp(-b2B*(RAB - c2B)*(RAB - c2B));
				}
				if (
					((coreMethod != "RM1-1g") && (coreMethod != "qRM1-1g"))
					&& (coreMethod[0] == 'q')
					)
				{
					F_A += a2A*expq((-b2A*(RAB - c2A)*(RAB - c2A)),qgauss);
					F_B += a2B*expq((-b2B*(RAB - c2B)*(RAB - c2B)),qgauss);
				}
				if ((coreMethod == "RM1-3g") || (coreMethod == "RM1"))
				{
					F_A += a3A*exp(-b3A*(RAB - c3A)*(RAB - c3A));
					F_B += a3B*exp(-b3B*(RAB - c3B)*(RAB - c3B));
				}
				if (coreMethod == "qRM1-3g")
				{
					F_A += a3A*expq((-b3A*(RAB - c3A)*(RAB - c3A)),qgauss);
					F_B += a3B*expq((-b3B*(RAB - c3B)*(RAB - c3B)),qgauss);
				}
			}
	
			za_zb = Params::get_int(mol.atom_name[A],"atomic_charge")*
				Params::get_int(mol.atom_name[B], "atomic_charge");
			gamass = Params::hartree_ev*four_center_.get_four_center(A, B, 0, 0, 0, 0);

			double expalfa;
			bool isExchange = (mol.scfMethod == "xRHF");
			double qalfa = Params::get_double(mol.atom_name[A], "qalfa");
			if (
				(coreMethod[0] == 'q') &&
				(abs(qalfa - 1.0e0) > 0.000001) &&
				(!isExchange)
				)
			{

				if((qalfa < 1) && 
				  (  (alfA*RAB >= 1/(1-qalfa)) || (alfB*RAB >= 1/(1-qalfa))  ))
					expalfa = 0.0e0;
				else
					expalfa = expq((-alfA*RAB), qalfa) + expq((-alfB*RAB), qalfa);

			}	
			else
			{

				expalfa = exp(-alfA*RAB) + exp(-alfB*RAB);
			}
			core_core_repulsion += za_zb*gamass*(1.0e0 + expalfa);
			core_core_repulsion += za_zb*(F_A + F_B)/RAB;
		}
	}
        if( core_core_repulsion < 0.0e0 )
	{
		return 0.0e0;
	}
	else
		return core_core_repulsion;
}

double ScfCycle::expq(double x,double q)
{
	if (abs(q - 1.0e0) > 1.0e0)
	{
		double aux = pow((1.0e0 + (1.0e0 - q)*x), (1.0e0 / (1.0e0 - q)));
		if (isnan(aux))
		{
			initializationSucess = false;
			return 0.0e0;
		}
		else
			return aux;
	}
	else
	{
		return exp(x);
	}
}

double ScfCycle::gaussSkew(double r,double a, double b, double c, double d)
{
	double pi = 4.0e0 * atan(1);
	double xdis = (d * ((r-c)/b))/sqrt(2.0e0);
	double erfunc = tanh(sqrt(pi) * log(2.0e0) * xdis);
	double skew = 0.5e0 * (1.0e0 + erfunc);
	double xgauss = (r-c)/b;
	double gauss = (1/sqrt(2.0e0 * pi)) * exp(-xgauss * xgauss/2.0e0);
	double gSkew  = a * (2.0e0/b) * gauss  * skew;
	return gSkew;
}

void ScfCycle::scfCycleInitialization()
{
	four_center_.build_matrix_integrals(mol);
	if (!four_center_.notANumber)
	{
		pPrintLog_->printNonZeroIntegrals(mol.number_of_atoms, four_center_.matrix_integrals, mol.atom_name);
		CoreMatrixCalculations coreMatrix_(mol, four_center_);
		if(coreMatrix_.getCoreSucess())
		{
			coreFockMatrix = coreMatrix_.fock_matrix;
			matrDiag_.symmetrizeEntryMatrix(coreFockMatrix);
			pPrintLog_->printScfMatrix("core", coreFockMatrix, Params::hartree_ev);
			fock_matrix_size = coreFockMatrix.size();
			num_eletrons = mol.number_of_electrons;
			scfMeth_.startScfMethod(scfMethod, fock_matrix_size, mol, four_center_, pPrintLog_);
			initializationSucess = true;
		}
		else
		{
			initializationSucess = false;
		}
	}
	else
	{
		initializationSucess = false;
	}
}

void ScfCycle::printScfIteration(int iStep)
{
	double elecEnergy = scfMeth_.getElectronicEnergy();
	if (iStep != 1)
	{
		double dE = elecEnergy - oldElectronicEnergy;
		double rmsDensity = calculateRmsDensity();
		pPrintLog_->printIteration(iStep, 
			elecEnergy,
			elecEnergy-oldElectronicEnergy,
			rmsDensity);
		
		if ((dE < 1.6e-7) && (rmsDensity < 1.6e-7))
		{
			if (checkDiis())
				converged = true;
		}
	}
	else
	{
		pPrintLog_->printScfHeader();
		pPrintLog_->printIteration(iStep, elecEnergy);
		oldDensity = scfMeth_.getDensity();
	}
	oldElectronicEnergy = elecEnergy;
}

double ScfCycle::calculateRmsDensity()
{
	vector< vector<double> > actualDensity = scfMeth_.getDensity();
	int size = actualDensity.size();
	double auxSum = 0.0e0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			auxSum += sqrt((actualDensity[i][j] - oldDensity[i][j])
				*(actualDensity[i][j] - oldDensity[i][j]));
		}
	}
	oldDensity = actualDensity;
	return auxSum/((double)size*size);
}

bool ScfCycle::checkDiis()
{
	if (scfMeth_.getDiis())
	{
		scfMeth_.setDiis(false);
		pPrintLog_->printEndDiis();
		return false;
	}
	else
	{
		return true;
	}
}

double ScfCycle::getIonizationPotential()
{
	return scfMeth_.getIonizationPotential();
}



