#ifndef MATRIXROTATION_H
#define MATRIXROTATION_H

#include<math.h>
#include <vector>

#include "Molecule.h"

class MatrixRotation
{
public:
	std::vector< std::vector<double> > CalcRotatingMatrix(const Molecule& mol, int A, int B);
	std::vector< std::vector<double> > CalcRotatingMatrix_over(const Molecule& mol, int A, int B);

};


#endif