#include "OptimizeWithDlib.h"

#include <dlib/optimization.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "FunctionDlib.h"
#include "FunctionDlibH2.h"
#include "FunctionDlibH3p.h"
#include "FunctionDlibn.h"
#include "PrintAll.h"
#include "Molecule.h"
#include "ScfProcedure.h"
#include "Params.h"
#include "ZmatoCart.h"

using namespace std;
using namespace dlib;

OptimizeWithDlib::~OptimizeWithDlib(){}

void OptimizeWithDlib::optimizen(
	PrintAll &printLog_, 
	Molecule &mol, 
	double finalEnergy,
	std::vector<int> &conections,
	std::vector<double> &startingPoint,
	ScfProcedure &scf_
	)
{
	int zmatSize = startingPoint.size();
	int size = (6 + zmatSize) / 3;
	column_vector starting_point(zmatSize);
	for (int i = 0; i < zmatSize; i++)
	{
		starting_point(i) = startingPoint[i];
	}

	find_min_using_approximate_derivatives(
		bfgs_search_strategy(),
		objective_delta_stop_strategy(1.0e-6
		),
		FunctionDlibn(&mol, printLog_, conections),
		starting_point,
		finalEnergy
		);

	// ENERGIA LIMITE ADMINISTRAR 50%

	//repetir o calculo.
	std::vector<CoordZMAT> zmat;
	zmat.resize(size);
	zmat[1].dis = starting_point(0);
	int k = 0;
	if (size > 2)
	{
		zmat[2].dis = starting_point(1);
		zmat[2].ang = starting_point(2);
		zmat[2].conect[0] = conections[k]; k++;
		zmat[2].conect[1] = conections[k]; k++;
	}
	int l = 2;
	for (int i = 3; i < size; i++)
	{
		l++; zmat[i].dis = starting_point(l);
		l++; zmat[i].ang = starting_point(l);
		l++; zmat[i].die = starting_point(l);
		zmat[i].conect[0] = conections[k]; k++;
		zmat[i].conect[1] = conections[k]; k++;
		zmat[i].conect[2] = conections[k]; k++;
	}
	ZmatoCart zmt_;	
	zmt_.ztocart(zmat);
	printLog_.printString("  ULTIMO  ");
	scf_.startScfProcedure(&mol, &printLog_);
	scf_.doScfProcedure(zmt_.MolXYZ);

	//quero retornar esse scf

	/*
	scf_.startScfProcedure(&mol, &printLog_);
	xyz[1] = starting_point(0);
	printLog_.printString("  ULTIMO  ");
	double totalEnergy = scf_.doScfProcedure(xyz);
	double ionizPot = scf_.getIonizationPotential();

	double error1 = abs(Params::ev_hartree*totalEnergy + 1.1727563206e0) / 1.1727563206e0;
	double error2 = abs(ionizPot - 0.6025790895e0) / 0.6025790895e0;
	double error3 = abs(xyz[1] - 0.7431918e0) / 0.7431918e0;

	double error = (error1 + error2 + error3) / 3.0e0;

	ofstream errorFile;
	errorFile.open("fitness-h2.txt");
	errorFile << error;
	errorFile.close();

	if (printExcell)
	{
		ofstream excell_;
		excell_.open("excell.txt");
		excell_ << setprecision(8) << Params::ev_hartree*totalEnergy << endl;
		excell_ << setprecision(8) << ionizPot << endl;
		excell_ << setprecision(8) << xyz[1] << endl;
		excell_.close();
	}
	*/
}

