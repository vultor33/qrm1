#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <dlib/optimization.h>

#include "Molecule.h"
#include "FourCenterIntegrals.h"
#include "CoreMatrixCalculations.h"
#include "ScfCycle.h"
#include "MatrixDiagonalization.h"
#include "ZmatoCart.h"
#include "Params.h"
#include "MatrixRotation.h"
#include "ZmatoCart.h"
#include "DiatomicOverlaps.h"
#include "PopleOverlaps.h"
#include "ScfProcedure.h"
#include "OptimizeWithDlib.h"
#include "ReadqInput.h"
#include "PrintAll.h"
#include "FitnessOpt.h"

using namespace std;
using namespace dlib;


/*
ATENCAO - mudanas para funcionar o q
- fredmudar em:
double FourCenterIntegrals::do_the_calculation_of_the_integral(const Molecule &mol, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma)

*/


int main(int argc, char *argv[])
{
	string qinputName;
	if (argc == 1)
		qinputName = "qinput.txt";
	else
		qinputName = argv[1];

	ReadqInput readQ_(qinputName.c_str());
	Params::set_semiempirical_parameters("RM1");
	Params::readFromFileParametrization("hparam.txt");
	PrintAll printLog_(readQ_.printLogName() , readQ_.getDebugLevel());
	printLog_.printOpening();
	printLog_.prinParameters();
	Molecule mol(readQ_.getCharge(), readQ_.getLabels());
	mol.setPrintLog(printLog_);
	mol.scfMethod = readQ_.getScfType();
	ScfProcedure scf_;


	/* fredmudar - testes pontuais
	scf_.startScfProcedure(&mol, &printLog_);
	std::vector<double> coord0 = readQ_.getCoordinates();
	scf_.doScfProcedure(coord0);
	cout << "PARAMETERIZATION DESACTIVATED CHECK SOURCE CODE" << endl;
	exit(1);
	*/


	OptimizeWithDlib opt_;
	std::vector<int> conections;
	std::vector<double> startingPoint;
	double finalEnergy;
	FitnessOpt fit_;
	switch (readQ_.getParametricType())
	{
	case 1:
		printLog_.printString("Modelo: H2p");
		conections.resize(0);
		startingPoint.resize(1);
		startingPoint[0] = 0.992155993;
		finalEnergy = -24.58591895;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(),"h2p");
		break;

	case 2:
		printLog_.printString("Modelo: H2");
		conections.resize(0);
		startingPoint.resize(1);
		startingPoint[0] = 0.697846959;
		finalEnergy = -47.8669996;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h2");
		break;

	case 3:
		printLog_.printString("Modelo: H3p");
		conections.resize(2);
		conections[0] = 1;
		conections[1] = 2;
		startingPoint.resize(3);
		startingPoint[0] = 0.868342;
		startingPoint[1] = 0.82093;
		startingPoint[2] = 62.46;
		finalEnergy = -54.7714018;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h3p");
		break;

	case 4:
		printLog_.printString("Modelo: H4p");
		conections.resize(5);
		conections[0] = 1;
		conections[1] = 2;
		conections[2] = 3;
		conections[3] = 1;
		conections[4] = 2;
		startingPoint.resize(6);
		startingPoint[0] = 0.773728575;
		startingPoint[1] = 0.8309298;
		startingPoint[2] = 58.32814636;
		startingPoint[3] = 1.872897344;
		startingPoint[4] = 147.3263107;
		startingPoint[5] = 167.58;
		finalEnergy = -75.48172424;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h4p");
		break;

	case 5:
		printLog_.printString("Modelo: H5p");
		conections.resize(8);
		conections[0] = 1;
		conections[1] = 2;
		conections[2] = 3;
		conections[3] = 1;
		conections[4] = 2;
		conections[5] = 4;
		conections[6] = 3;
		conections[7] = 1;
		startingPoint.resize(9);
		startingPoint[0] = 0.779337364;
		startingPoint[1] = 0.9911196;
		startingPoint[2] = 66.04598;
		startingPoint[3] = 1.186160096;
		startingPoint[4] = 150.5497142;
		startingPoint[5] = 129.023709;
		startingPoint[6] = 0.7047145;
		startingPoint[7] = 66.47712031;
		startingPoint[8] = 123.1866081;
		finalEnergy = -103.1716149;
		opt_.optimizen(printLog_, mol, finalEnergy,conections,startingPoint,scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h5p");
		break;

	case 6:
		printLog_.printString("Modelo: H6p");
		conections.resize(11);
		conections[0] = 2;
		conections[1] = 1;
		conections[2] = 3;
		conections[3] = 2;
		conections[4] = 1;
		conections[5] = 1;
		conections[6] = 2;
		conections[7] = 3;
		conections[8] = 5;
		conections[9] = 1;
		conections[10] = 2;
		startingPoint.resize(12);
		startingPoint[0] = 1.131735;
		startingPoint[1] = 1.103125315;
		startingPoint[2] = 153.7998121;
		startingPoint[3] = 0.86193;
		startingPoint[4] = 64.59267255;
		startingPoint[5] = 163.98;
		startingPoint[6] = 1.271336689;
		startingPoint[7] = 157.3465583;
		startingPoint[8] = -91.62;
		startingPoint[9] = 0.75831;
		startingPoint[10] = 66.72914463;
		startingPoint[11] = 174.6;
		finalEnergy = -123.9943096;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h6p");
		break;

	case 7:
		printLog_.printString("Modelo: H7p");
		conections.resize(14);
		conections[0] = 1;
		conections[1] = 2;
		conections[2] = 3;
		conections[3] = 1;
		conections[4] = 2;
		conections[5] = 4;
		conections[6] = 3;
		conections[7] = 1;
		conections[8] = 2;
		conections[9] = 1;
		conections[10] = 3;
		conections[11] = 6;
		conections[12] = 2;
		conections[13] = 1;
		startingPoint.resize(15);
		startingPoint[0] = 0.8217066;
		startingPoint[1] = 0.8400249;
		startingPoint[2] = 64.9312;
		startingPoint[3] = 1.689713814;
		startingPoint[4] = 147.1175911;
		startingPoint[5] = 159.3708941;
		startingPoint[6] = 0.7038552;
		startingPoint[7] = 72.5218353;
		startingPoint[8] = 114.7916321;
		startingPoint[9] = 1.541185095;
		startingPoint[10] = 149.2475252;
		startingPoint[11] = -148.5965238;
		startingPoint[12] = 0.7974008;
		startingPoint[13] = 71.14773737;
		startingPoint[14] = -107.9839494;
		finalEnergy = -153.748486;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h7p");
		break;

	case 8:
		printLog_.printString("Modelo: H8p");
		conections.resize(17);
		conections[0] = 2;
		conections[1] = 1;
		conections[2] = 3;
		conections[3] = 2;
		conections[4] = 1;
		conections[5] = 1;
		conections[6] = 2;
		conections[7] = 3;
		conections[8] = 5;
		conections[9] = 1;
		conections[10] = 2;
		conections[11] = 3;
		conections[12] = 2;
		conections[13] = 1;
		conections[14] = 7;
		conections[15] = 3;
		conections[16] = 2;
		startingPoint.resize(18);
		startingPoint[0] = 0.95304;
		startingPoint[1] = 1.114672;
		startingPoint[2] = 154.8396835;
		startingPoint[3] = 0.767845875;
		startingPoint[4] = 72.41098048;
		startingPoint[5] = -164.5199991;
		startingPoint[6] = 1.316812302;
		startingPoint[7] = 148.5960972;
		startingPoint[8] = -66.8462889;
		startingPoint[9] = 0.724768;
		startingPoint[10] = 72.01952634;
		startingPoint[11] = 155.7394194;
		startingPoint[12] = 2.282928932;
		startingPoint[13] = 121.0170215;
		startingPoint[14] = 11.29374602;
		startingPoint[15] = 0.784719;
		startingPoint[16] = 81.02468108;
		startingPoint[17] = -102.2782833;
		finalEnergy = -171.9820808;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h8p");
		break;

	case 9:
		printLog_.printString("Modelo: H9p");
		conections.resize(20);
		conections[0] = 1;
		conections[1] = 2;
		conections[2] = 3;
		conections[3] = 1;
		conections[4] = 2;
		conections[5] = 4;
		conections[6] = 3;
		conections[7] = 1;
		conections[8] = 2;
		conections[9] = 1;
		conections[10] = 3;
		conections[11] = 6;
		conections[12] = 2;
		conections[13] = 1;
		conections[14] = 1;
		conections[15] = 2;
		conections[16] = 3;
		conections[17] = 8;
		conections[18] = 1;
		conections[19] = 2;
		startingPoint.resize(21);
		startingPoint[0] = 0.9052317;
		startingPoint[1] = 0.9746595;
		startingPoint[2] = 61.98;
		startingPoint[3] = 1.667624463;
		startingPoint[4] = 140.4021612;
		startingPoint[5] = 167.3224622;
		startingPoint[6] = 0.6761104;
		startingPoint[7] = 82.88808901;
		startingPoint[8] = 100.3838751;
		startingPoint[9] = 1.713522384;
		startingPoint[10] = 147.3410683;
		startingPoint[11] = -146.3098748;
		startingPoint[12] = 0.8156848;
		startingPoint[13] = 75.93568639;
		startingPoint[14] = -107.5937658;
		startingPoint[15] = 1.718622153;
		startingPoint[16] = 152.803614;
		startingPoint[17] = 158.9174268;
		startingPoint[18] = 0.7819168;
		startingPoint[19] = 77.78966043;
		startingPoint[20] = 121.5698643;
		finalEnergy = -199.5014923;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h9p");
		break;

	default:
		printLog_.printString("METODO DE OTIMIZACAO NAO ENCONTRADO");
		break;
	}

}


