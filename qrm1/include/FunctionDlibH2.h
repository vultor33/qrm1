#ifndef FUNCTIONDLIBH2_H
#define FUNCTIONDLIBH2_H

#include <dlib/optimization.h>
#include <vector>

#include "Molecule.h"
#include "PrintAll.h"
#include "ScfProcedure.h"
#include "ZmatoCart.h"
#include "Coordstructs.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlibH2
{
public:
	FunctionDlibH2(Molecule * pMol_in, PrintAll &printLog_in_)
	{
		scf_.startScfProcedure(pMol_in, &printLog_in_);
		size = pMol_in->number_of_atoms;
	}

	~FunctionDlibH2(){}

	// h2
	double operator() (const column_vector& arg) const
	{
		std::vector<double> xyz(6);
		for (int i = 0; i < 6; i++)
		{
			xyz[i] = 0.0e0;
		}
		xyz[1] = arg(0);
		ScfProcedure scf2_ = scf_;
		return scf2_.doScfProcedure(xyz);
	}

private:
	ScfProcedure scf_;
	ZmatoCart zmt_;
	int size;


};

#endif