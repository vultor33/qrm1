#ifndef FUNCTIONDLIB_H
#define FUNCTIONDLIB_H

#include <dlib/optimization.h>
#include <vector>

#include "Molecule.h"
#include "PrintAll.h"
#include "ScfProcedure.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlib
{
public:
	FunctionDlib(Molecule * pMol_in, PrintAll &printLog_in_)
	{
		scf_.startScfProcedure(pMol_in, &printLog_in_);
	}

	~FunctionDlib(){}

	double operator() (const column_vector& arg) const
	{
		std::vector<double> xyz(arg.size());
		for (int i = 0; i < arg.size(); i++)
		{
			xyz[i] = arg(i);
		}
		ScfProcedure scf2_ = scf_;
		return scf2_.doScfProcedure(xyz);
	}

private:
	ScfProcedure scf_;


};

#endif