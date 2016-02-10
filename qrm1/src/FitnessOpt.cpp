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
	errorFile << error;
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
		refEnergy = -0.602344;
		refIonizPot = 1.102636;
		refDistances.resize(1);
		refDistances[0] = 1.0577356;
		break;

	case 2:
		refEnergy = -1.172719;
		refIonizPot = 0.602545;
		refDistances.resize(1);
		refDistances[0] = 0.743181;
		break;

	case 3:
		refEnergy = -1.341874;
		refIonizPot = 1.226814;
		refDistances.resize(3);
		refDistances[0] = 0.878000;
		refDistances[1] = 0.878000;
		refDistances[2] = 0.878000;
		break;

	case 4:
		refEnergy = -1.849267;
		refIonizPot = 0.748213;
		refDistances.resize(6);
		refDistances[0] = 0.845605;
		refDistances[1] = 0.890600;
		refDistances[2] = 0.890600;
		refDistances[3] = 2.556785;
		refDistances[4] = 2.522567;
		refDistances[5] = 1.721413;
		break;

	case 5:
		refEnergy = -2.527657;
		refIonizPot = 0.9031266368;
		refDistances.resize(10);
		refDistances[0] = 0.8101219;
		refDistances[1] = 0.9951;
		refDistances[2] = 0.9951;
		refDistances[3] = 2.234407;
		refDistances[4] = 2.234407;
		refDistances[5] = 1.3121239;
		refDistances[6] = 2.234407;
		refDistances[7] = 2.234407;
		refDistances[8] = 1.312124;
		refDistances[9] = 0.768500;
		break;

	case 6:
		refEnergy = -3.037803;
		refIonizPot = 0.706506;
		refDistances.resize(15);
		refDistances[0] = 1.045000;
		refDistances[1] = 2.233754;
		refDistances[2] = 1.218923;
		refDistances[3] = 2.233754;
		refDistances[4] = 1.218923;
		refDistances[5] = 0.785000;
		refDistances[6] = 1.2189226;
		refDistances[7] = 2.2337541;
		refDistances[8] = 3.3986352;
		refDistances[9] = 3.3986352;
		refDistances[10] = 1.2189226;
		refDistances[11] = 2.2337541;
		refDistances[12] = 3.3986352;
		refDistances[13] = 3.3986352;
		refDistances[14] = 0.785;
		break;

	case 7:
		refEnergy = -3.766767;
		refIonizPot = 0.924027;
		refDistances.resize(21);
		refDistances[0] = 0.872300;
		refDistances[1] = 0.872300;
		refDistances[2] = 0.949156;
		refDistances[3] = 2.403094;
		refDistances[4] = 2.436317;
		refDistances[5] = 1.597083;
		refDistances[6] = 2.4030935;
		refDistances[7] = 2.4363167;
		refDistances[8] = 1.5970828;
		refDistances[9] = 0.7544;
		refDistances[10] = 2.4030935;
		refDistances[11] = 1.5970828;
		refDistances[12] = 2.4363167;
		refDistances[13] = 3.5662916;
		refDistances[14] = 3.6452099;
		refDistances[15] = 2.4030935;
		refDistances[16] = 1.5970828;
		refDistances[17] = 2.4363167;
		refDistances[18] = 3.6452099;
		refDistances[19] = 3.5662916;
		refDistances[20] = 0.7544;
		break;

	case 8:
		refEnergy = -4.213482;
		refIonizPot = 0.695968;
		refDistances.resize(28);
		refDistances[0] = 1.045000;
		refDistances[1] = 2.210265;
		refDistances[2] = 1.196000;
		refDistances[3] = 2.180026;
		refDistances[4] = 1.167000;
		refDistances[5] = 0.781125;
		refDistances[6] = 1.2505336;
		refDistances[7] = 2.2637071;
		refDistances[8] = 3.3894606;
		refDistances[9] = 3.3885971;
		refDistances[10] = 1.2505336;
		refDistances[11] = 2.2637071;
		refDistances[12] = 3.3894606;
		refDistances[13] = 3.3885971;
		refDistances[14] = 0.781;
		refDistances[15] = 3.5172305;
		refDistances[16] = 2.8667014;
		refDistances[17] = 2.1060225;
		refDistances[18] = 2.8740452;
		refDistances[19] = 4.3532686;
		refDistances[20] = 4.4191853;
		refDistances[21] = 3.4990354;
		refDistances[22] = 2.8551933;
		refDistances[23] = 2.1055095;
		refDistances[24] = 2.8744616;
		refDistances[25] = 4.3966611;
		refDistances[26] = 4.3304016;
		refDistances[27] = 0.741;
		break;

	case 9:
		refEnergy = -4.887695;
		refIonizPot = 0.806144;
		refDistances.resize(36);
		refDistances[0] = 0.890100;
		refDistances[1] = 0.890100;
		refDistances[2] = 0.890100;
		refDistances[3] = 2.497630;
		refDistances[4] = 2.497630;
		refDistances[5] = 1.699923;
		refDistances[6] = 2.4976295;
		refDistances[7] = 2.4976295;
		refDistances[8] = 1.6999232;
		refDistances[9] = 0.7504;
		refDistances[10] = 2.4976295;
		refDistances[11] = 1.6999232;
		refDistances[12] = 2.4976295;
		refDistances[13] = 3.7618402;
		refDistances[14] = 3.8359539;
		refDistances[15] = 2.4976295;
		refDistances[16] = 1.6999232;
		refDistances[17] = 2.4976295;
		refDistances[18] = 3.8359539;
		refDistances[19] = 3.7618402;
		refDistances[20] = 0.7504;
		refDistances[21] = 1.6999232;
		refDistances[22] = 2.4976295;
		refDistances[23] = 2.4976295;
		refDistances[24] = 3.7618402;
		refDistances[25] = 3.8359539;
		refDistances[26] = 3.7618402;
		refDistances[27] = 3.8359539;
		refDistances[28] = 1.6999232;
		refDistances[29] = 2.4976295;
		refDistances[30] = 2.4976295;
		refDistances[31] = 3.8359539;
		refDistances[32] = 3.7618402;
		refDistances[33] = 3.8359539;
		refDistances[34] = 3.7618402;
		refDistances[35] = 0.7504;
		break;

	default:
		cout << "Erro em FitnessOpt::setSystems - modelo n encontrado" << endl;
		break;
	}
}

