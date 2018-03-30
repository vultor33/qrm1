#include "PrintAll.h"

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Coordstructs.h"
#include "Molecule.h"
#include "FourCenterMatrix.h"
#include "Params.h"
#include "Coordstructs.h"

using namespace std;

PrintAll::PrintAll(string outputName_in, int debugLevel)
{
	outputName = outputName_in;
	string logName = outputName + ".log";
	logOutput_.open(logName.c_str());

	debugOptions(debugLevel);
}

PrintAll::~PrintAll()
{
	logOutput_.close();
}

void PrintAll::debugOptions(int debugLevel)
{
	if (debugLevel == 0)
	{
		conditionFirst = false;
		conditionMatrixIterations = false;
		conditionXyzCoordinates = true;
		conditionDiis = false;
	}
	else
	{
		conditionFirst = true;
		conditionMatrixIterations = true;
		conditionXyzCoordinates = true;
		conditionDiis = true;
	}
}

bool PrintAll::checkPrintConditions(string type)
{
	if ((type == "core") && (conditionFirst))
		return true;
	if ((type == "first density") && (conditionFirst))
		return true;
	if ((type == "non zero integrals") && (conditionFirst))
		return true;
	if ((type == "fock") && (conditionMatrixIterations))
		return true;
	if ((type == "density") && (conditionMatrixIterations))
		return true;
	if ((type == "diis error") && (conditionDiis))
		return true;
	if ((type == "diis b coefficients") && (conditionDiis))
		return true;
	if ((type == "diis old fock") && (conditionDiis))
		return true;
	if ((type == "diis new fock") && (conditionDiis))
		return true;
	if ((type == "diis error") && (conditionDiis))
		return true;
	if ((type == "xyz") && (conditionXyzCoordinates))
		return true;

	return false;
}

void PrintAll::printString(string anyString)
{
	logOutput_ << anyString << endl;
}

void PrintAll::printOpening()
{
	logOutput_ << "================================" << endl
		       << "    q-semiempirical package     " << endl
		       << "================================" << endl
		       << "              author:Fred Vultor" << endl << endl;
}

void PrintAll::prinParameters()
{
	logOutput_ << "METHOD:  " << Params::method << endl
		<< "H parameters: " << endl
		<< "Uss:  " << Params::get_double("H", "uss") << endl
		<< "Betas:  " << Params::get_double("H", "betas") << endl
		<< "alfacore:  " << Params::get_double("H", "alfacore") << endl
		<< "gss:  " << Params::get_double("H", "gss") << endl
		<< "a1core:  " << Params::get_double("H", "a1core") << endl
		<< "b1core:  " << Params::get_double("H", "b1core") << endl
		<< "c1core:  " << Params::get_double("H", "c1core") << endl
		<< "zetas:  " << Params::get_double("H", "expoents") << endl
		<< "qover:  " << Params::get_double("H", "qover") << endl
		<< "qmono:  " << Params::get_double("H", "qmono") << endl
		<< "qalfa:  " << Params::get_double("H", "qalfa") << endl
		<< "qgauss:  " << Params::get_double("H", "qgauss") << endl;
	if((Params::method == "RM1-2g") || (Params::method == "qRM1-2g"))
	{
		logOutput_ << "a2core:  " << Params::get_double("H", "a2core") << endl
		<< "b2core:  " << Params::get_double("H", "b2core") << endl
		<< "c2core:  " << Params::get_double("H", "c2core") << endl;
	}
	else if((Params::method == "RM1-3g") || (Params::method == "qRM1-3g"))
	{
		logOutput_ << "a2core:  " << Params::get_double("H", "a2core") << endl
		<< "b2core:  " << Params::get_double("H", "b2core") << endl
		<< "c2core:  " << Params::get_double("H", "c2core") << endl
		<< "a3core:  " << Params::get_double("H", "a3core") << endl
		<< "b3core:  " << Params::get_double("H", "b3core") << endl
		<< "c3core:  " << Params::get_double("H", "c3core") << endl;
	}
	logOutput_ << endl << endl;

}

void PrintAll::printMatrix(const vector< vector<double> > &entryMatrix, double unitConversion)
{
	int lines = entryMatrix.size();
	int cols = entryMatrix[0].size();
	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			logOutput_ << setiosflags(ios::fixed) << setprecision(6) << setw(12) <<
				entryMatrix[i][j]*unitConversion << "  ";
		}
		logOutput_ << endl;
	}
	logOutput_ << endl;
}


void PrintAll::printScfMatrix(
	string matrixType,
	const vector< vector<double> > &entryMatrix,
	double unitConversion
	)
{
	if (checkPrintConditions(matrixType))
	{
		logOutput_ << endl;
		logOutput_ << matrixType << endl;
		printMatrix(entryMatrix, unitConversion);
	}
}

void PrintAll::printScfMatrix(
	string matrixType,
	const vector<double> &entryMatrix,
	double unitConversion
	)
{
	if (checkPrintConditions(matrixType))
	{
		logOutput_ << matrixType << endl;
		int size = entryMatrix.size();
		int jumpLine = 0;
		for (int i = 0; i < size; i++)
		{
			logOutput_ << setw(8) << setprecision(6)
				<< entryMatrix[i] << "    ";
			jumpLine++;
			if (jumpLine == 5)
			{
				jumpLine = 0;
				logOutput_ << endl;
			}
		}
		logOutput_ << endl << endl;
	}
}

void PrintAll::printMolecule(
	const vector<CoordZMAT> & zmatMolecule,
	const vector<CoordXYZ> & xyzMolecule,
	int nAtoms
	)
{
	// debug sobre o que imprimir e tal
//	logOutput_ << "INPUT MOLECULE " << endl;
//	logOutput_ << "Coordinates (mopac style) " << endl;
//	printZmatOnLogOutput(zmatMolecule, nAtoms);

	if (checkPrintConditions("xyz"))
	{
		logOutput_ << "INPUT MOLECULE " << endl;
		logOutput_ << "Coordinates (angs) " << endl;
		printXyzOnLogOutput(xyzMolecule, nAtoms);
	}
}

void PrintAll::printOnXyzFile(const vector<CoordXYZ> & xyzMolecule, int nAtoms)
{
	ofstream xyzFile_;
	string xyzName = outputName + ".xyz";
	xyzFile_.open(xyzName.c_str());

	xyzFile_ << nAtoms << endl << endl;

	for (int i = 0; i<nAtoms; i++)
	{
		xyzFile_ << xyzMolecule[i].atomlabel << "  "
			<< xyzMolecule[i].x << "   "
			<< xyzMolecule[i].y << "   "
			<< xyzMolecule[i].z << endl;
	}
}

void PrintAll::printXyzOnLogOutput(const vector<CoordXYZ> & xyzMolecule, int nAtoms)
{
	logOutput_ << "           X         Y          Z  " << endl;
	for (int i = 0; i<nAtoms; i++)
	{
		logOutput_ 
			<< xyzMolecule[i].atomlabel << "  "
			<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
			<< xyzMolecule[i].x << "   "
			<< xyzMolecule[i].y << "   "
			<< xyzMolecule[i].z << endl;
	}
	logOutput_ << endl;
}

void PrintAll::printZmatOnLogOutput(const vector<CoordZMAT> & zmatMolecule, int nAtoms)
{
	for (int i = 0; i < nAtoms; i++)
	{
		logOutput_ << setiosflags(ios::fixed) << setw(4)
			<< zmatMolecule[i].atomlabel << "   "
			<< setprecision(6)
			<< zmatMolecule[i].dis << "   "
			<< setw(3)
			<< zmatMolecule[i].flagconec[0] << "   "
			<< setprecision(6)
			<< zmatMolecule[i].ang << "   "
			<< setw(3)
			<< zmatMolecule[i].flagconec[1] << "   "
			<< setprecision(6)
			<< zmatMolecule[i].die << "   "
			<< setw(3)
			<< zmatMolecule[i].flagconec[2] << "   ";
		if (i > 0)
			logOutput_ << setw(3) << zmatMolecule[i].conect[0] << "   ";
		if (i > 1)
			logOutput_ << setw(3) << zmatMolecule[i].conect[1] << "   ";
		if (i > 2)
			logOutput_ << setw(3) << zmatMolecule[i].conect[2];
		
		logOutput_ << endl;
	}
	logOutput_ << endl;
}


void PrintAll::printNonZeroIntegrals(
	int nAtoms,
	vector< vector<diatomicFourCenter> > &matrixIntegrals,
	vector< string > & atomName
	)
{
	if (checkPrintConditions("non zero integrals"))
	{
		logOutput_ << "NON ZERO INTEGRALS (eV)" << endl;

		int fullLine = 0;
		for (int A = 0; A < nAtoms; A++)
		{
			int nBasesA = Params::get_int(atomName[A], "base_number");
			for (int B = A; B < nAtoms; B++)
			{
				int nBasesB = Params::get_int(atomName[B], "base_number");
				for (int mi = 0; mi < nBasesA; mi++)
				{
					for (int ni = 0; ni < nBasesB; ni++)
					{
						for (int lambda = 0; lambda < nBasesA; lambda++)
						{
							for (int sigma = 0; sigma < nBasesB; sigma++)
							{
								if (abs(matrixIntegrals[A][B].f_center[mi][ni][lambda][sigma]) > 1.0e-6)
								{
									logOutput_ << setiosflags(ios::fixed) << setprecision(6) << setw(12)
										<< matrixIntegrals[A][B].f_center[mi][ni][lambda][sigma] * Params::hartree_ev
										<< "    ";
									fullLine++;
									if (fullLine == 4)
									{
										logOutput_ << endl;
										fullLine = 0;
									}
								}
							}
						}
					}
				}
			}
		}
		logOutput_ << endl << endl;
	}
}

void PrintAll::printIteration(int iStep, double energy, double energyVariation, double densityRms)
{
	logOutput_ 
		<< setiosflags(ios::fixed) << setprecision(6)
		<< setw(4) << iStep 
		<< "     " << setw(10) << energy*Params::hartree_ev
		<< "     " << setw(10) << energyVariation*Params::hartree_ev
		<< "     " << setw(10) << densityRms
		<< endl;
}

void PrintAll::printScfHeader(int flag)
{
	if ((flag == 0)||conditionDiis)
	{
		logOutput_ << "SELF CONSISTENT FIELD" << endl
			<< "Iter       E                dE             Rms(D)" << endl;
	}
}

void PrintAll::printFinalEnergy(
	double elecEnergy, 
	double coreEnergy, 
	double ionizationPotential)
{
	logOutput_  << "Total energy:  " 
		<< fixed << setprecision(8) << coreEnergy + elecEnergy << endl
		<< "Final electronic energy:  " 
		<< fixed << setprecision(8) << elecEnergy << endl
		<< "Core repulsion energy:  " 
		<< fixed << setprecision(8) << coreEnergy << endl
		<< "Ionization potential:  "
		<< fixed << setprecision(8) << ionizationPotential << endl
		<< endl << endl;
}

void PrintAll::printEndOfScf(bool converged)
{
	if (converged)
	{
		logOutput_ << endl
			<< " SCF CONVERGED SUCCESSFULLY" << endl << endl;
	}
	else
	{
		logOutput_ << endl
			<< " SCF DIDN'T CONVERGED" << endl << endl;
	}
}

void PrintAll::printStartDiis()
{
	logOutput_ << "    DIIS PROCEDURE TURNED ON       " << endl;
}

void PrintAll::printErrorDiis(double maxError)
{
	if(checkPrintConditions("diis error"))
		logOutput_ << "Max diss error matrix:  " << maxError << endl;
}


void PrintAll::printEndDiis()
{
	logOutput_ << "    DIIS PROCEDURE TURNED OFF       " << endl;
}

