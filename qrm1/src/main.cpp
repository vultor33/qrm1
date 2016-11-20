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

int main(int argc, char *argv[])
{
	string qinputName;
	if (argc == 1)
		qinputName = "qinput.txt";
	else
		qinputName = argv[1];

	ReadqInput readQ_(qinputName.c_str());
	Params::set_semiempirical_parameters("RM1");
	PrintAll printLog_(readQ_.printLogName() , readQ_.getDebugLevel());
	printLog_.printOpening();
	Molecule mol(readQ_.getCharge(), readQ_.getLabels());
	mol.setPrintLog(printLog_);
	mol.scfMethod = readQ_.getScfType();
	Params::readFromFileParametrization("hparam.txt");

	OptimizeWithDlib opt_;
	std::vector<int> conections;
	std::vector<double> startingPoint;
	double finalEnergy;
	ScfProcedure scf_;

	//fredmudar
	scf_.doScfProcedure("h3.xyz");


	FitnessOpt fit_;
	switch (readQ_.getParametricType())
	{
	case 1:
		printLog_.printString("Modelo: H2p");
		conections.resize(0);
		startingPoint.resize(1);
		startingPoint[0] = 1.0577356;
		finalEnergy = -24.58591895;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(),"h2p");
		break;

	case 2:
		printLog_.printString("Modelo: H2");
		conections.resize(0);
		startingPoint.resize(1);
		startingPoint[0] = 0.743181;
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
		startingPoint[0] = 0.878000;
		startingPoint[1] = 0.878000;
		startingPoint[2] = 60.000000;
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
		startingPoint[0] = 0.845605;
		startingPoint[1] = 0.890600;
		startingPoint[2] = 61.65766;
		startingPoint[3] = 1.721413;
		startingPoint[4] = 155.080327;
		startingPoint[5] = 180.000000;
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
		startingPoint[0] = 0.810122;
		startingPoint[1] = 0.995100;
		startingPoint[2] = 65.980000;
		startingPoint[3] = 1.312124;
		startingPoint[4] = 150.851417;
		startingPoint[5] = 143.041806;
		startingPoint[6] = 0.768500;
		startingPoint[7] = 72.971592;
		startingPoint[8] = 123.309918;
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
		startingPoint[0] = 1.045000;
		startingPoint[1] = 1.218923;
		startingPoint[2] = 161.215736;
		startingPoint[3] = 0.785000;
		startingPoint[4] = 71.215736;
		startingPoint[5] = 180.000000;
		startingPoint[6] = 1.218923;
		startingPoint[7] = 161.215736;
		startingPoint[8] = -90.000000;
		startingPoint[9] = 0.785000;
		startingPoint[10] = 71.215736;
		startingPoint[11] = 180.000000;
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
		startingPoint[0] = 0.872300;
		startingPoint[1] = 0.872300;
		startingPoint[2] = 65.920000;
		startingPoint[3] = 1.597083;
		startingPoint[4] = 152.138150;
		startingPoint[5] = 149.644032;
		startingPoint[6] = 0.754400;
		startingPoint[7] = 76.338774;
		startingPoint[8] = 117.373857;
		startingPoint[9] = 1.597083;
		startingPoint[10] = 152.138150;
		startingPoint[11] = -149.644032;
		startingPoint[12] = 0.754400;
		startingPoint[13] = 76.338774;
		startingPoint[14] = -117.373858;
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
		startingPoint[0] = 1.045000;
		startingPoint[1] = 1.196000;
		startingPoint[2] = 160.956012;
		startingPoint[3] = 0.781125;
		startingPoint[4] = 68.701120;
		startingPoint[5] = -179.999999;
		startingPoint[6] = 1.250534;
		startingPoint[7] = 160.818287;
		startingPoint[8] = -71.877730;
		startingPoint[9] = 0.781000;
		startingPoint[10] = 71.804114;
		startingPoint[11] = 160.887830;
		startingPoint[12] = 2.106023;
		startingPoint[13] = 117.835464;
		startingPoint[14] = 11.465732;
		startingPoint[15] = 0.741000;
		startingPoint[16] = 79.827272;
		startingPoint[17] = -93.064862;
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
		startingPoint[0] = 0.890100;
		startingPoint[1] = 0.890100;
		startingPoint[2] = 60.000000;
		startingPoint[3] = 1.699923;
		startingPoint[4] = 147.636342;
		startingPoint[5] = 155.648802;
		startingPoint[6] = 0.750400;
		startingPoint[7] = 77.248918;
		startingPoint[8] = 110.921409;
		startingPoint[9] = 1.699923;
		startingPoint[10] = 147.636341;
		startingPoint[11] = -155.648803;
		startingPoint[12] = 0.750400;
		startingPoint[13] = 77.248918;
		startingPoint[14] = -110.921408;
		startingPoint[15] = 1.699923;
		startingPoint[16] = 147.636342;
		startingPoint[17] = 155.648802;
		startingPoint[18] = 0.750400;
		startingPoint[19] = 77.248918;
		startingPoint[20] = 110.921409;
		finalEnergy = -199.5014923;
		opt_.optimizen(printLog_, mol, finalEnergy, conections, startingPoint, scf_);
		fit_.printFitness(scf_, readQ_.getPrintExcell(), readQ_.getParametricType(), "h9p");
		break;

	default:
		printLog_.printString("METODO DE OTIMIZACAO NAO ENCONTRADO");
		break;
	}

}


