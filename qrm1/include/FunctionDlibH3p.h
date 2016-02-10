#ifndef FUNCTIONDLIBH3P_H
#define FUNCTIONDLIBH3P_H

#include <dlib/optimization.h>
#include <vector>

#include "Molecule.h"
#include "PrintAll.h"
#include "ScfProcedure.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlibH3p
{
public:
	FunctionDlibH3p(Molecule * pMol_in, PrintAll &printLog_in_)
	{
		scf_.startScfProcedure(pMol_in, &printLog_in_);
	}

	~FunctionDlibH3p(){}

	// h3
	double operator() (const column_vector& arg) const
	{
		std::vector<double> xyz(9);
		for (int i = 0; i < 9; i++)
		{
			xyz[i] = 0.0e0;
		}
		xyz[1] = arg(0);
		double auxy = arg(1)*sin(arg(2));
		double auxxA = xyz[1];
		double auxx = arg(1)*cos(arg(2)) + auxxA;
		xyz[2] = auxx;
		xyz[5] = auxy;

		ScfProcedure scf2_ = scf_;
		return scf2_.doScfProcedure(xyz);
	}

private:
	ScfProcedure scf_;
};

#endif

