#include "ScfMethod.h"

#include <iostream>

#include "RhfMethod.h"
#include "UhfMethod.h"

using namespace std;

ScfMethod::~ScfMethod()
{
	if (sucessfulStart)
	{
		if (
			//exchangemudar
			(method == "RHF")||
			(method == "xRHF")
			)
		{
			rhfMatrix_.deletePointers();
		}
		else
		{
			uhfMatrix_.deletePointers();
		}
	}
}

void ScfMethod::startScfMethod(string method_in, int fock_matrix_size, Molecule & mol, FourCenterMatrix & four_center_, PrintAll * pPrintLog_in_)
{
	method = method_in;

	if (
		//exchangemudar
		(method == "RHF")||
		(method == "xRHF")
		)
	{
		rhfMatrix_.startRhfFockMatrix(fock_matrix_size, mol, four_center_, pPrintLog_in_);
		rhfMatrix_.build_first_density_matrix();
	}
	else
	{
		uhfMatrix_.startUhfFockMatrix(fock_matrix_size, mol, four_center_, pPrintLog_in_);
		uhfMatrix_.build_first_density_matrix();
	}
	sucessfulStart = true;
}

void ScfMethod::scfMethodCycle(int diagMethod, vector< vector<double> > &coreFockMatrix)
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		rhfMatrix_.setDiagonalizationMethod(diagMethod);
		rhfMatrix_.buildFockMatrix(coreFockMatrix);
		rhfMatrix_.electronic_energy(coreFockMatrix);
	}
	else
	{
		uhfMatrix_.setDiagonalizationMethod(diagMethod);
		uhfMatrix_.buildFockMatrix(coreFockMatrix);
		uhfMatrix_.electronic_energy(coreFockMatrix);
	}
}

double ScfMethod::getElectronicEnergy()
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		return rhfMatrix_.getElectronicEnergy();
	}
	else
	{
		return uhfMatrix_.getElectronicEnergy();
	}
}

vector< vector<double> > ScfMethod::getDensity()
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		return rhfMatrix_.getDensity();
	}
	else
	{
		return uhfMatrix_.getDensity();
	}
}

bool ScfMethod::getDiis()
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		return rhfMatrix_.getDiis();
	}
	else
	{
		return uhfMatrix_.getDiis();
	}
}


void ScfMethod::setDiis(bool newDiis)
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		return rhfMatrix_.setDiis(newDiis);
	}
	else
	{
		return uhfMatrix_.setDiis(newDiis);
	}
}



double ScfMethod::getIonizationPotential()
{
	if (
		//exchangemudar
		(method == "RHF") ||
		(method == "xRHF")
		)
	{
		return rhfMatrix_.getIonizationPotential();
	}
	else
	{
		return uhfMatrix_.getIonizationPotential();
	}
}










