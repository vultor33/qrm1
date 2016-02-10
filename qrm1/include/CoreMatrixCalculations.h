#ifndef COREMATRIXCALCULATIONS_H
#define COREMATRIXCALCULATIONS_H

#include <vector>

#include "FourCenterIntegrals.h"
#include "FourCenterMatrix.h"
#include "DiatomicOverlaps.h"

class CoreMatrixCalculations
{
private:
	double mi_mi_integral(const Molecule &mol, int A, int mi);
	double mi_ni_integral(const Molecule &mol, int A, int mi, int ni);
	double mi_ni_potential(const Molecule &mol, int A, int B, int mi, int ni);
	void set_fock_matrix_dimension(Molecule &mol);
	void print_core_matrix();
        bool notANumber;

	FourCenterMatrix &four_center_;
	FourCenterIntegrals four_center_integrals_;
	DiatomicOverlaps ressonance_AB_;

public:
	CoreMatrixCalculations(Molecule &mol, FourCenterMatrix &four_center_in);
	void build_fock_matrix(Molecule &mol);
	std::vector< std::vector<double> > fock_matrix;

	inline bool getCoreSucess(){ return !notANumber;}

};

#endif
