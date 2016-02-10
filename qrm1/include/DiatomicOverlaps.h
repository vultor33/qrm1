#ifndef DIATOMICOVERLAPS_H
#define DIATOMICOVERLAPS_H

#include <vector>

#include "Molecule.h"
#include "CheckRessonanceType.h"
#include "Matrixop.h"
#include "MatrixRotation.h"
#include "PopleOverlaps.h"

struct diatomic_overlap_AB
{
	std::vector< std::vector<double> > bases;
};

class DiatomicOverlaps
{
	std::vector< std::vector<diatomic_overlap_AB> > overlaps_AB;

	std::vector< std::vector<double> > calc_overlap_AB(const Molecule &mol, int A, int B);
	void overlapToRessonance(const Molecule &mol, std::vector< std::vector<double> > &over_AB, int A, int B);
	double getBetaAtom(const Molecule &mol, int orbital, int Atom);
	void pople_rotation(const Molecule &mol, std::vector< std::vector<double> > &over_AB, int A, int B);
	double calc_overlap_AB_mini(const Molecule &mol, int A, int B, int mi, int ni);

	double overlap_1s_1s(const Molecule &mol, int A, int B);
	double overlap_1s_2s(const Molecule &mol, int A, int B);
	double overlap_1s_2pz(const Molecule &mol, int A, int B);
	double overlap_2s_2s(const Molecule &mol, int A, int B);
	double overlap_2s_2pz(const Molecule &mol, int A, int B);
	double overlap_2pz_2pz(const Molecule &mol, int A, int B);
	double overlap_2ppi_2ppi(const Molecule &mol, int A, int B);

	CheckRessonanceType check_;
	Matrixop matrix_operation_;
	MatrixRotation m_rotation_;
	PopleOverlaps popOver_;


public:
//	DiatomicOverlaps(const Molecule &mol);
	void start_overlaps(const Molecule &mol);

	double get_overlap_AB(int A, int B, int mi, int ni);

};



#endif