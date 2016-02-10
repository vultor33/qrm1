#include "ReadqInput.h"

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

ReadqInput::~ReadqInput(){}

ReadqInput::ReadqInput(string inputName)
{
	ifstream qinp_;
	string dummy;
	qinp_.open(inputName.c_str());

	qinp_ >> dummy >> parametricType;
	//cout << "parametricType:  " << parametricType << endl;

	qinp_ >> dummy >> nAtoms;
	//cout << "natomos:  " << nAtoms << endl;

	qinp_ >> dummy >> charge;
	//cout << "charge:  " << charge << endl;

	qinp_ >> dummy >> scfType;
	//cout << "scfType:  " << scfType << endl;

	qinp_ >> dummy >> debugLevel;
	//cout << "debugLevel:  " << debugLevel << endl;

	int auxPrintExcell;
	qinp_ >> dummy >> auxPrintExcell;
	//	cout << "printExcell:  " << printExcell << endl;
	if (auxPrintExcell == 0)
		printExcell = false;
	else
		printExcell = true;

	labels.resize(nAtoms);
	coordinates.resize(nAtoms * 3);
	for (int i = 0; i < nAtoms; i++)
	{
		qinp_ >> labels[i];
//		cout << "labels[i]:  " << labels[i] << endl;
	}
	for (int i = 0; i < 3*nAtoms; i++)
	{
		qinp_ >> coordinates[i];
//		cout << "coordinates[i]:  " << coordinates[i] << endl;
	}

	qinp_.close();
}

std::string ReadqInput::printLogName()
{
	switch (parametricType)
	{
	case 1:
		return "rm1-h2p";
	case 2:
		return "rm1-h2";
	case 3:
		return "rm1-h3p";
	case 4:
		return "rm1-h4p";
	case 5:
		return "rm1-h5p";
	case 6:
		return "rm1-h6p";
	case 7:
		return "rm1-h7p";
	case 8:
		return "rm1-h8p";
	case 9:
		return "rm1-h9p";
	default:
		return "ERRO-PRINTLOGNAME";
	}
}

std::string ReadqInput::getScfType()
{
	switch (scfType)
	{
	case 0:
		return "RHF";
	case 1:
		return "UHF";
	case 2:
		return "xRHF";
	default:
		cout << "ERRO EM getScfTYPE tipo de scf nao cadastrado" << endl;
		exit(3);
	}
}








