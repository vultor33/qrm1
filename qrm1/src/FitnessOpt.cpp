#include "FitnessOpt.h"

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "ScfProcedure.h"
#include "Params.h"
#include "Molecule.h"

using namespace std;

FitnessOpt::~FitnessOpt(){}

void FitnessOpt::printFitness(ScfProcedure &scf_, bool printExcell, int model, string errorName)
{
	setSystems(model);
	double totalEnergy = scf_.getFinalEnergy();
	double ionizPot = scf_.getIonizationPotential();
	Molecule *pMol = scf_.getPmol();

	double error1 = abs(Params::ev_hartree*totalEnergy - refEnergy) / abs(refEnergy);
	double error2 = abs(ionizPot - refIonizPot) / refIonizPot;
	double error3 = calculateDistanceError(pMol);

	double error = (error1 + error2 + error3) / 3.0e0;

	ofstream errorFile;
	errorName += ".ga";
	errorFile.open(errorName.c_str());
	errorFile << setprecision(16) << error;
	errorFile.close();

	if (printExcell)
	{
		ofstream excell_;
		excell_.open("excell.txt");
		excell_ << setprecision(8) << Params::ev_hartree*totalEnergy << endl;
		excell_ << setprecision(8) << ionizPot << endl;
		printExcellDistances(excell_, pMol);
		excell_.close();
	}
}

double FitnessOpt::calculateDistanceError(Molecule *pMol)
{
	int size = pMol->number_of_atoms;
	int k = 0;
	double error = 0.0e0;
	for(int j = 1; j < size; j++)
	{
		for (int i = 0; i < j; i++)
		{
			error += abs(pMol->distancia_i_j[i][j] * Params::bohr_angs - refDistances[k]) / refDistances[k];
			k++;
		}
	}
	error /= k;
	return error;
}

void FitnessOpt::printExcellDistances(ofstream &excell_, Molecule *pMol)
{
	int size = pMol->number_of_atoms;
	for (int j = 1; j < size; j++)
	{
		for (int i = 0; i < j; i++)
		{
			excell_ << setprecision(8) << pMol->distancia_i_j[i][j] * Params::bohr_angs << endl;
		}
	}
}

void FitnessOpt::setSystems(int model)
{
	switch (model)
	{
	case 1:
		refEnergy = -0.602344e0;
		refIonizPot = 1.102636e0;
		refDistances.resize(1);
		refDistances[0] = 1.0577356e0;
		break;

	case 2:
		refEnergy = -1.172719e0;
		refIonizPot = 0.602545e0;
		refDistances.resize(1);
		refDistances[0] = 0.743181e0;
		break;

	case 3:
		refEnergy = -1.341874e0;
		refIonizPot = 1.226814e0;
		refDistances.resize(3);
		refDistances[0] = 0.878000e0;
		refDistances[1] = 0.878000e0;
		refDistances[2] = 0.878000e0;
		break;

	case 4:
		refEnergy = -1.849267e0;
		refIonizPot = 0.748213e0;
		refDistances.resize(6);
		refDistances[0] = 0.845605e0;
		refDistances[1] = 0.890600e0;
		refDistances[2] = 0.890600e0;
		refDistances[3] = 2.556785e0;
		refDistances[4] = 2.522567e0;
		refDistances[5] = 1.721413e0;
		break;

	case 5:
		refEnergy = -2.527657e0;
		refIonizPot = 0.9031266368e0;
		refDistances.resize(10);
		refDistances[0] = 0.8101219e0;
		refDistances[1] = 0.9951e0;
		refDistances[2] = 0.9951e0;
		refDistances[3] = 2.234407e0;
		refDistances[4] = 2.234407e0;
		refDistances[5] = 1.3121239e0;
		refDistances[6] = 2.234407e0;
		refDistances[7] = 2.234407e0;
		refDistances[8] = 1.312124e0;
		refDistances[9] = 0.768500e0;
		break;

	case 6:
		refEnergy = -3.037803e0;
		refIonizPot = 0.706506e0;
		refDistances.resize(15);
		refDistances[0] = 1.045000e0;
		refDistances[1] = 2.233754e0;
		refDistances[2] = 1.218923e0;
		refDistances[3] = 2.233754e0;
		refDistances[4] = 1.218923e0;
		refDistances[5] = 0.785000e0;
		refDistances[6] = 1.2189226e0;
		refDistances[7] = 2.2337541e0;
		refDistances[8] = 3.3986352e0;
		refDistances[9] = 3.3986352e0;
		refDistances[10] = 1.2189226e0;
		refDistances[11] = 2.2337541e0;
		refDistances[12] = 3.3986352e0;
		refDistances[13] = 3.3986352e0;
		refDistances[14] = 0.785e0;
		break;

	case 7:
		refEnergy = -3.766767e0;
		refIonizPot = 0.924027e0;
		refDistances.resize(21);
		refDistances[0] = 0.872300e0;
		refDistances[1] = 0.872300e0;
		refDistances[2] = 0.949156e0;
		refDistances[3] = 2.403094e0;
		refDistances[4] = 2.436317e0;
		refDistances[5] = 1.597083e0;
		refDistances[6] = 2.4030935e0;
		refDistances[7] = 2.4363167e0;
		refDistances[8] = 1.5970828e0;
		refDistances[9] = 0.7544e0;
		refDistances[10] = 2.4030935e0;
		refDistances[11] = 1.5970828e0;
		refDistances[12] = 2.4363167e0;
		refDistances[13] = 3.5662916e0;
		refDistances[14] = 3.6452099e0;
		refDistances[15] = 2.4030935e0;
		refDistances[16] = 1.5970828e0;
		refDistances[17] = 2.4363167e0;
		refDistances[18] = 3.6452099e0;
		refDistances[19] = 3.5662916e0;
		refDistances[20] = 0.7544e0;
		break;

	case 8:
		refEnergy = -4.213482e0;
		refIonizPot = 0.695968e0;
		refDistances.resize(28);
		refDistances[0] = 1.045000e0;
		refDistances[1] = 2.210265e0;
		refDistances[2] = 1.196000e0;
		refDistances[3] = 2.180026e0;
		refDistances[4] = 1.167000e0;
		refDistances[5] = 0.781125e0;
		refDistances[6] = 1.2505336e0;
		refDistances[7] = 2.2637071e0;
		refDistances[8] = 3.3894606e0;
		refDistances[9] = 3.3885971e0;
		refDistances[10] = 1.2505336e0;
		refDistances[11] = 2.2637071e0;
		refDistances[12] = 3.3894606e0;
		refDistances[13] = 3.3885971e0;
		refDistances[14] = 0.781e0;
		refDistances[15] = 3.5172305e0;
		refDistances[16] = 2.8667014e0;
		refDistances[17] = 2.1060225e0;
		refDistances[18] = 2.8740452e0;
		refDistances[19] = 4.3532686e0;
		refDistances[20] = 4.4191853e0;
		refDistances[21] = 3.4990354e0;
		refDistances[22] = 2.8551933e0;
		refDistances[23] = 2.1055095e0;
		refDistances[24] = 2.8744616e0;
		refDistances[25] = 4.3966611e0;
		refDistances[26] = 4.3304016e0;
		refDistances[27] = 0.741e0;
		break;

	case 9:
		refEnergy = -4.887695e0;
		refIonizPot = 0.806144e0;
		refDistances.resize(36);
		refDistances[0] = 0.890100e0;
		refDistances[1] = 0.890100e0;
		refDistances[2] = 0.890100e0;
		refDistances[3] = 2.497630e0;
		refDistances[4] = 2.497630e0;
		refDistances[5] = 1.699923e0;
		refDistances[6] = 2.4976295e0;
		refDistances[7] = 2.4976295e0;
		refDistances[8] = 1.6999232e0;
		refDistances[9] = 0.7504e0;
		refDistances[10] = 2.4976295e0;
		refDistances[11] = 1.6999232e0;
		refDistances[12] = 2.4976295e0;
		refDistances[13] = 3.7618402e0;
		refDistances[14] = 3.8359539e0;
		refDistances[15] = 2.4976295e0;
		refDistances[16] = 1.6999232e0;
		refDistances[17] = 2.4976295e0;
		refDistances[18] = 3.8359539e0;
		refDistances[19] = 3.7618402e0;
		refDistances[20] = 0.7504e0;
		refDistances[21] = 1.6999232e0;
		refDistances[22] = 2.4976295e0;
		refDistances[23] = 2.4976295e0;
		refDistances[24] = 3.7618402e0;
		refDistances[25] = 3.8359539e0;
		refDistances[26] = 3.7618402e0;
		refDistances[27] = 3.8359539e0;
		refDistances[28] = 1.6999232e0;
		refDistances[29] = 2.4976295e0;
		refDistances[30] = 2.4976295e0;
		refDistances[31] = 3.8359539e0;
		refDistances[32] = 3.7618402e0;
		refDistances[33] = 3.8359539e0;
		refDistances[34] = 3.7618402e0;
		refDistances[35] = 0.7504e0;
		break;

	default:
		cout << "Erro em FitnessOpt::setSystems - modelo n encontrado" << endl;
		break;
	}
}

