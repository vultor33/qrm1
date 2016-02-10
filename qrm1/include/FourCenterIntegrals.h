#ifndef FOURCENTERINTEGRALS_H
#define FOURCENTERINTEGRALS_H

#include <vector>

#include "Molecule.h"
#include "CheckTypeOfIntegralToBeCalculated.h"

class FourCenterIntegrals
{
private:
	double do_the_calculation_of_the_integral(const Molecule &mol, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);

	double q_q_bracket(const Molecule &mol, int atom_A, int atom_B);
	double q_miz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double q_Qzz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double q_Qpipi_bracket(const Molecule &mol, int atom_A, int atom_B);
	double mipi_mipi_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qpipi_Qpipi_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qxx_Qyy_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qpipi_Qzz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qzz_Qzz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double miz_Qpipi_bracket(const Molecule &mol, int atom_A, int atom_B);
	double miz_Qzz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double miz_miz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double mipi_Qpiz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qpiz_Qpiz_bracket(const Molecule &mol, int atom_A, int atom_B);
	double Qxy_Qxy_bracket(const Molecule &mol, int atom_A, int atom_B);

	double overlap_sign(int atom_A, int atom_B);

	// vou cadastrar todas as integrais assim, acho q e mais pratico de testar
	// ai depois eu tiro essa logistica.
	int check_if_it_is_calculated(int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	std::vector < std::vector<int> > flag_of_calculated_integrals;
	std::vector<double> calculated_integrals;

	CheckTypeOfIntegralToBeCalculated check_;

public:
	double get_four_center_integral(const Molecule &mol, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);

};

#endif