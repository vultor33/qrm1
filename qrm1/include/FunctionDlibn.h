#ifndef FUNCTIONDLIBN_H
#define FUNCTIONDLIBN_H

#include <dlib/optimization.h>
#include <vector>

#include "Molecule.h"
#include "PrintAll.h"
#include "ScfProcedure.h"
#include "ZmatoCart.h"
#include "Coordstructs.h"

typedef dlib::matrix<double, 0, 1> column_vector;

class FunctionDlibn
{
public:
	FunctionDlibn(Molecule * pMol_in, PrintAll &printLog_in_, std::vector<int> &conections_in)
	{
		scf_.startScfProcedure(pMol_in, &printLog_in_);
		size = pMol_in->number_of_atoms;
		conections = conections_in;
	}

	~FunctionDlibn(){}

	//conections: atom 3.1, 3.2; atom 4.1, 4.2, 4.3 e assim por diante
	//arg: dis 2; dis3, ang3; dis4, ang4, die4 e assim por diante
	double operator() (const column_vector& arg) const
	{
		std::vector<CoordZMAT> zmat;
		zmat.resize(size);
		zmat[1].dis = arg(0);
		int k = 0;
		if (size > 2)
		{
			zmat[2].dis = arg(1);
			zmat[2].ang = arg(2);
			zmat[2].conect[0] = conections[k]; k++;
			zmat[2].conect[1] = conections[k]; k++;
		}
		int l = 2;
		for (int i = 3 ; i < size; i++)
		{
			l++; zmat[i].dis = arg(l);
			l++; zmat[i].ang = arg(l);
			l++; zmat[i].die = arg(l);
			zmat[i].conect[0] = conections[k]; k++;
			zmat[i].conect[1] = conections[k]; k++;
			zmat[i].conect[2] = conections[k]; k++;
		}
		ZmatoCart zmt_;
		zmt_.ztocart(zmat);
		ScfProcedure scf2_ = scf_;
		return scf2_.doScfProcedure(zmt_.MolXYZ);
	}

private:
	ScfProcedure scf_;
	int size;
	std::vector<int> conections;

};

#endif





