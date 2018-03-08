#include "ScfProcedure.h"

#include <string>

#include "PrintAll.h"
#include "Molecule.h"
#include "Params.h"
#include "ScfCycle.h"

using namespace std;

ScfProcedure::~ScfProcedure(){}


void ScfProcedure::startScfProcedure(Molecule * mol_in, PrintAll * printLog_in_)
{
	pMol = mol_in;
	pPrintLog_ = printLog_in_;
}


void ScfProcedure::doScfProcedure(string inputName)
{
	// object that prints all output
	PrintAll printLog_("qRM1-teste", 1);
	printLog_.printOpening();

	// parameters setting
	Params::set_semiempirical_parameters("RM1");
//	Params::readFromFile("hparam.txt");

	// input reading
	Molecule mol(inputName);
	mol.setPrintLog(printLog_);

	// internal to cartesian and distances
	mol.build_atoms();

	// self consistent field
	string scfMethod;
	if (mol.number_of_electrons % 2 == 0)
		scfMethod = "RHF";
	else
		scfMethod = "UHF";

	ScfCycle cycle_(mol, printLog_, scfMethod);

	if (cycle_.qNotNan())
	{
		bool converged;
		for (int i = 1; i < 50; i++)
		{
			converged = cycle_.SCF_step_i(i);
			if (converged)
				break;
		}
		printLog_.printEndOfScf(converged);
		cycle_.getFinalEnergy();
		cout << "Work Complete" << endl;
	}
	else
	{
		cout << "Din't start well" << endl;
	}
}

//double ScfProcedure::doScfProcedure(const column_vector& arg)
double ScfProcedure::doScfProcedure(vector<double> &xyzCoord)
{
	// setting cooordinates
	pMol->buildXyzAtoms(xyzCoord);

	// self consistent field
	ScfCycle cycle_(*pMol, *pPrintLog_, pMol->scfMethod);
	if (cycle_.qNotNan())
	{
		bool converged;
		for (int i = 1; i < 100; i++)
		{
			converged = cycle_.SCF_step_i(i);
			if (converged||!(cycle_.qNotNan()))
				break;
		}
		if (cycle_.qNotNan())
		{
			pPrintLog_->printEndOfScf(converged);
			ionizationPotential = cycle_.getIonizationPotential();
			if (converged)
			{
				finalEnergy = cycle_.getFinalEnergy();
				return finalEnergy;
			}
			else
			{
				finalEnergy = 0.0e0;
				return finalEnergy;
			}
		}
	}
	finalEnergy = 1.0e0;
	return finalEnergy;
}



//double ScfProcedure::doScfProcedure(const column_vector& arg)
double ScfProcedure::doScfProcedure(vector<CoordXYZ> &xyzCoord)
{
	// setting cooordinates
	pMol->buildXyzAtoms(xyzCoord);

	// self consistent field
	ScfCycle cycle_(*pMol, *pPrintLog_, pMol->scfMethod);
	if (cycle_.qNotNan())
	{
		bool converged;
		for (int i = 1; i < 50; i++)
		{
			converged = cycle_.SCF_step_i(i);
			if (converged || !(cycle_.qNotNan()))
				break;
		}
		if (cycle_.qNotNan())
		{
			pPrintLog_->printEndOfScf(converged);
			ionizationPotential = cycle_.getIonizationPotential();
			if (converged)
			{
				finalEnergy = cycle_.getFinalEnergy();
				return finalEnergy;
			}
			else
			{
				finalEnergy = 0.0e0;
				return finalEnergy;
			}
		}
	}
	finalEnergy = 1.0e0;
	return finalEnergy;
}

