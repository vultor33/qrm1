#ifndef FOURCENTERMATRIX_H
#define FOURCENTERMATRIX_H

#include <vector>

#include "FourCenterIntegrals.h"
#include "Molecule.h"
#include "MatrixRotation.h"
#include "Coordstructs.h"

class FourCenterMatrix
{
public:
	void build_matrix_integrals(Molecule &mol);

	bool notANumber;
	std::vector< std::vector<diatomicFourCenter> > matrix_integrals;
	double get_four_center(int A, int B, int mi, int ni, int lambda, int sigma, bool isExchange = false);


private:
	void four_center_rotation(const Molecule &mol, int A,int B);
	int transforma(int base);

	Molecule * pMol_;
	FourCenterIntegrals integrals_;
	MatrixRotation m_rotation_;
};

#endif

