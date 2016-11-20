#include "FourCenterIntegrals.h"

#include<iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "Molecule.h"
#include "Params.h"
#include "CheckTypeOfIntegralToBeCalculated.h"

using namespace std;

double FourCenterIntegrals::get_four_center_integral(const Molecule &mol, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma)
{
	int calculou = check_if_it_is_calculated(atomo_A,atomo_B,orbital_mi,orbital_ni,orbital_lambda,orbital_sigma);
	if (calculou == 0)
	{
		return do_the_calculation_of_the_integral(mol,atomo_A,atomo_B,orbital_mi,orbital_ni,orbital_lambda,orbital_sigma);
	}
	else
	{
		return calculated_integrals[calculou];
	}
}

int FourCenterIntegrals::check_if_it_is_calculated(int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma)
{
	int int_calc = flag_of_calculated_integrals.size();
	if (int_calc == 0)
	{
		return 0;
	}
	else
	{
		for (int i = 0; i < int_calc; i++)
		{
			int A = flag_of_calculated_integrals[i][0];
			int B = flag_of_calculated_integrals[i][1];
			int mi = flag_of_calculated_integrals[i][2];
			int ni = flag_of_calculated_integrals[i][3];
			int lambda = flag_of_calculated_integrals[i][4];
			int sigma = flag_of_calculated_integrals[i][5];

			// Condicoes para que a integral tenha sido calculada.
			// existem mais, mas em primeiro momento apenas
			// essas serão consideradas.
			if (

				(
				(A == atomo_A) && (B == atomo_B)
				) &&

				(
				((mi == orbital_mi) && (ni == orbital_ni))
				||
				((mi == orbital_ni) && (ni == orbital_mi))
				) &&

				(
				((lambda == orbital_lambda) && (sigma == orbital_sigma))
				||
				((lambda == orbital_sigma) && (sigma == orbital_lambda))
				)

				)
			{
				return i;
			}
		}
	}
	return 0;
}


double FourCenterIntegrals::do_the_calculation_of_the_integral(const Molecule &mol, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma)
{
	double integral_AB;
//	cout << "four_center::do_the_calculation  A:  " << atomo_A << "  B:  " << atomo_B << "  mi:  " << orbital_mi << "  ni  " << orbital_ni << " lam  " << orbital_lambda << "  sig  " << orbital_sigma << endl;
	if (atomo_A == atomo_B)
	{
		switch (check_.check_one_center_type(orbital_mi, orbital_ni, orbital_lambda, orbital_sigma))
		{
		case 0:
			return 0.0e0;
		case 1:
			//cout << "ssss A:  " << atomo_A << "  B:  " << atomo_B << "  mi:  " << orbital_mi << "  ni  " << orbital_ni << " lam  " << orbital_lambda << "  sig  " << orbital_sigma << endl;
			//fredmudar
			//return Params::get_double(mol.atom_name[atomo_A], "qr0");
			return Params::get_double(mol.atom_name[atomo_A], "gss");
		case 2:
			return Params::get_double(mol.atom_name[atomo_A], "gsp");
		case 3:
			return Params::get_double(mol.atom_name[atomo_A], "hsp");
		case 4:
			return Params::get_double(mol.atom_name[atomo_A], "gpp");
		case 5:
			return Params::get_double(mol.atom_name[atomo_A], "gp2");
		case 6:
			return 0.5e0*(Params::get_double(mol.atom_name[atomo_A], "gpp")
				- Params::get_double(mol.atom_name[atomo_A], "gp2"));
		default:
			cout << "nao sei calcular essa integral - erro em do_the_calculation_of_the_integral" << endl;
			exit(7);
		}
	}
	else
	{
		//cout << "four_center::do_the_calculation  A:  " << atomo_A << "  B:  " << atomo_B << "  mi:  " << orbital_mi << "  ni  " << orbital_ni << " lam  " << orbital_lambda << "  sig  " << orbital_sigma << endl;
		switch (check_.check_two_center_type(orbital_mi, orbital_ni, orbital_lambda, orbital_sigma))
		{
		case 1: //(ss|ss)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B);
			break;

		case 2: //(ss|spz) or (ss|pzs)
			integral_AB = q_miz_bracket(mol, atomo_A, atomo_B);
			break;

		case 3: //(spz|ss) or (pzs|ss)
			integral_AB = -q_miz_bracket(mol, atomo_B, atomo_A);
			break;

		case 4: //(ss|pzpz)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qzz_bracket(mol, atomo_A, atomo_B);
			break;

		case 5: //(pzpz|ss)
			integral_AB = q_q_bracket(mol, atomo_B, atomo_A) +
				q_Qzz_bracket(mol, atomo_B, atomo_A);
		break;

		case 6: //(ss|ppippi)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_A, atomo_B);
			break;

		case 7: //(ppippi|ss)
			integral_AB = q_q_bracket(mol, atomo_B, atomo_A) +
				q_Qpipi_bracket(mol, atomo_B, atomo_A);
			break;

		case 8: //(sppi|sppi)
			integral_AB = mipi_mipi_bracket(mol, atomo_A, atomo_B);
			break;

		case 9: //(ppippi|ppippi)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_B, atomo_A) +
				Qpipi_Qpipi_bracket(mol, atomo_A, atomo_B);
			break;

		case 10: //(pxpx|pypy)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_B, atomo_A) +
				Qxx_Qyy_bracket(mol, atomo_A, atomo_B);
			break;

		case 11: //(ppippi|pzpz)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qzz_bracket(mol, atomo_A, atomo_B) +
				q_Qpipi_bracket(mol, atomo_B, atomo_A) +
				Qpipi_Qzz_bracket(mol, atomo_A, atomo_B);
			break;

		case 12: //(pzpz|ppippi)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qzz_bracket(mol, atomo_B, atomo_A) +
				q_Qpipi_bracket(mol, atomo_A, atomo_B) +
				Qpipi_Qzz_bracket(mol, atomo_B, atomo_A);
			break;

		case 13: //(pzpz|pzpz)
			integral_AB = q_q_bracket(mol, atomo_A, atomo_B) +
				q_Qzz_bracket(mol, atomo_A, atomo_B) +
				q_Qzz_bracket(mol, atomo_B, atomo_A) +
				Qzz_Qzz_bracket(mol, atomo_A, atomo_B);
			break;

		case 14: //(spz|ppippi)
			integral_AB = - q_miz_bracket(mol, atomo_B, atomo_A) +
				miz_Qpipi_bracket(mol, atomo_A, atomo_B);
			break;

		case 15: //(ppippi|spz)
			integral_AB = q_miz_bracket(mol, atomo_A, atomo_B) -
				miz_Qpipi_bracket(mol, atomo_B, atomo_A);
			break;

		case 16: //(spz|pzpz)
			integral_AB = - q_miz_bracket(mol, atomo_B, atomo_A) +
				miz_Qzz_bracket(mol, atomo_A, atomo_B);
			break;

		case 17: //(pzpz|spz)
			integral_AB = q_miz_bracket(mol, atomo_A, atomo_B) -
				miz_Qzz_bracket(mol, atomo_B, atomo_A);
			break;

		case 18: //(spz|spz)
			integral_AB = miz_miz_bracket(mol, atomo_A, atomo_B);
			break;

		case 19: //(sppi|ppipz)
			integral_AB = mipi_Qpiz_bracket(mol, atomo_A, atomo_B);
			break;

		case 20: //(ppipz|sppi)
			integral_AB = - mipi_Qpiz_bracket(mol, atomo_B, atomo_A);
			break;

		case 21: //(ppipz|ppipz)
			integral_AB = Qpiz_Qpiz_bracket(mol, atomo_A, atomo_B);
			break;

		case 22: //(pxpy|pxpy)
			integral_AB = 0.5e0*(Qpipi_Qpipi_bracket(mol, atomo_A, atomo_B)-Qxx_Qyy_bracket(mol, atomo_A, atomo_B));
			break;

		default:
			return 0.0e0;
			//cout << "nao sei calcular essa integral - checar FourCenterIntegrals::do_the_calculation_of_the_integral" << endl;
			//exit(4);
		}
		//cout << "four center integral:  " << integral_AB*Params::hartree_ev << endl;

		return integral_AB;
	}
}

double FourCenterIntegrals::overlap_sign(int A,int B)
{
	if(A>B)
	{
		return -1.0e0;
	}
	else
	{
		return 1.0e0;
	}
}


double FourCenterIntegrals::q_q_bracket(const Molecule &mol, int atomo_A, int atomo_B)
{
	double R_distancia = mol.distancia_i_j[atomo_A][atomo_B];
	double ro_0_A = Params::get_double(mol.atom_name[atomo_A], "ro_0");
	double ro_0_B = 1.0e0 / (2.0e0 * Params::get_double(mol.atom_name[atomo_B], "gss"));
	double ro_0_B2 = Params::get_double(mol.atom_name[atomo_B], "ro_0");

	double q = Params::get_double(mol.atom_name[atomo_A], "qmono");
	double integral_AB;
	if (
		(Params::method[0] == 'q') &&
		(abs(q-1.0e0) > 0.000001)
		)
	{
		double x = R_distancia*R_distancia + (ro_0_A + ro_0_B)*(ro_0_A + ro_0_B);
		integral_AB = pow(1.5e0 - 0.5e0*pow(x, 1.0e0 - q), (1.0e0 / (1.0e0 - q)))
			- Params::get_double(mol.atom_name[atomo_A],"qrinf");
	}
	else
	{
		integral_AB = 1.0e0 / (sqrt(R_distancia*R_distancia + (ro_0_A + ro_0_B)*(ro_0_A + ro_0_B)));
	}

	return integral_AB;
}

double FourCenterIntegrals::q_miz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double ro_0_A = Params::get_double(mol.atom_name[A],"ro_0");
	double D1_B = Params::get_double(mol.atom_name[B], "D1");
	double ro_1_B = Params::get_double(mol.atom_name[B], "ro_1");

	double a_0_1_quad = (ro_0_A + ro_1_B)*(ro_0_A + ro_1_B);

	double colchete = 0.5e0*((1.0e0 / sqrt((RAB + D1_B)*(RAB + D1_B) + a_0_1_quad)) -
		1.0e0 / sqrt((RAB - D1_B)*(RAB - D1_B) + a_0_1_quad));

	return colchete;
}

double FourCenterIntegrals::q_Qzz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double ro_0_A = Params::get_double(mol.atom_name[A], "ro_0");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_0_2_quad = (ro_0_A + ro_2_B)*(ro_0_A + ro_2_B);

	double colchete = 0.25e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_B)*(RAB + 2.0e0*D2_B) +
		a_0_2_quad)) - 0.5e0*(1.0e0 / sqrt(RAB*RAB + a_0_2_quad)) +
		0.25e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_B)*(RAB - 2.0e0*D2_B) + a_0_2_quad));

	return colchete;
}


double FourCenterIntegrals::q_Qpipi_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double ro_0_A = Params::get_double(mol.atom_name[A], "ro_0");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_0_2_quad = (ro_0_A + ro_2_B)*(ro_0_A + ro_2_B);

	double colchete = 0.5e0*(1.0e0/sqrt(RAB*RAB+4.0e0*D2_B*D2_B+a_0_2_quad))
		- 0.5e0*(1.0e0 / sqrt(RAB*RAB+a_0_2_quad));

	return colchete;
}

double FourCenterIntegrals::mipi_mipi_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D1_A = Params::get_double(mol.atom_name[A], "D1");
	double ro_1_A = Params::get_double(mol.atom_name[A], "ro_1");
	double D1_B = Params::get_double(mol.atom_name[B], "D1");
	double ro_1_B = Params::get_double(mol.atom_name[B], "ro_1");

	double a_1_1_quad = (ro_1_A + ro_1_B)*(ro_1_A + ro_1_B);

	double colchete = 0.5e0*(
		(1.0e0 / sqrt(RAB*RAB + (D1_A - D1_B)*(D1_A - D1_B) + a_1_1_quad)) -
		(1.0e0 / sqrt(RAB*RAB + (D1_A + D1_B)*(D1_A + D1_B) + a_1_1_quad))
		);

	return colchete;

}

double FourCenterIntegrals::Qpipi_Qpipi_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = Params::get_double(mol.atom_name[A], "D2");
	double ro_2_A = Params::get_double(mol.atom_name[A], "ro_2");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete = 0.125e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*(D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		+ 0.125e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*(D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))
		- 0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_A*D2_A + a_2_2_quad))
		- 0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_B*D2_B + a_2_2_quad))
		+ 0.25e0*(1.0e0 / sqrt(RAB*RAB + a_2_2_quad));

	return colchete;
}

double FourCenterIntegrals::Qxx_Qyy_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = Params::get_double(mol.atom_name[A], "D2");
	double ro_2_A = Params::get_double(mol.atom_name[A], "ro_2");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete = 0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_A*D2_A + 4.0e0*D2_B*D2_B + a_2_2_quad))
		-0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_A*D2_A + a_2_2_quad))
		-0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_B*D2_B + a_2_2_quad))
		+0.25e0*(1.0e0 / sqrt(RAB*RAB + a_2_2_quad));

	return colchete;
}

double FourCenterIntegrals::Qpipi_Qzz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = Params::get_double(mol.atom_name[A], "D2");
	double ro_2_A = Params::get_double(mol.atom_name[A], "ro_2");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete = 0.125e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_B)*(RAB - 2.0e0*D2_B) + 4.0e0*D2_A*D2_A + a_2_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_B)*(RAB + 2.0e0*D2_B) + 4.0e0*D2_A*D2_A + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_B)*(RAB - 2.0e0*D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_B)*(RAB + 2.0e0*D2_B) + a_2_2_quad))
		- 0.25e0*(1.0e0 / sqrt(RAB*RAB + 4.0e0*D2_A*D2_A + a_2_2_quad))
		+ 0.25e0*(1.0e0 / sqrt(RAB*RAB + a_2_2_quad));

	return colchete;
}

double FourCenterIntegrals::Qzz_Qzz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = Params::get_double(mol.atom_name[A], "D2");
	double ro_2_A = Params::get_double(mol.atom_name[A], "ro_2");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete = 0.0625e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_A - 2.0e0*D2_B)*(RAB + 2.0e0*D2_A - 2.0e0*D2_B) + a_2_2_quad))
		+ 0.0625e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_A + 2.0e0*D2_B)*(RAB + 2.0e0*D2_A + 2.0e0*D2_B) + a_2_2_quad))
		+ 0.0625e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_A - 2.0e0*D2_B)*(RAB - 2.0e0*D2_A - 2.0e0*D2_B) + a_2_2_quad))
		+ 0.0625e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_A + 2.0e0*D2_B)*(RAB - 2.0e0*D2_A + 2.0e0*D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_A)*(RAB + 2.0e0*D2_A) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_A)*(RAB - 2.0e0*D2_A) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + 2.0e0*D2_B)*(RAB + 2.0e0*D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB - 2.0e0*D2_B)*(RAB - 2.0e0*D2_B) + a_2_2_quad))
		+ 0.25e0*(1.0e0 / sqrt(RAB*RAB + a_2_2_quad));

	return colchete;
}



double FourCenterIntegrals::miz_Qpipi_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D1_A = Params::get_double(mol.atom_name[A], "D1");
	double ro_1_A = Params::get_double(mol.atom_name[A], "ro_1");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_1_2_quad = (ro_1_A + ro_2_B)*(ro_1_A + ro_2_B);

	double colchete = - 0.25e0*(1.0e0 / sqrt((RAB + D1_A)*(RAB + D1_A) + 4.0e0*D2_B*D2_B + a_1_2_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB - D1_A)*(RAB - D1_A) + 4.0e0*D2_B*D2_B + a_1_2_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB + D1_A)*(RAB + D1_A) + a_1_2_quad))
		- 0.25e0*(1.0e0 / sqrt((RAB - D1_A)*(RAB - D1_A) + a_1_2_quad));

	return colchete;
}


double FourCenterIntegrals::miz_Qzz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D1_A = Params::get_double(mol.atom_name[A], "D1");
	double ro_1_A = Params::get_double(mol.atom_name[A], "ro_1");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_1_2_quad = (ro_1_A + ro_2_B)*(ro_1_A + ro_2_B);

	double colchete = -0.125e0*(1.0e0 / sqrt((RAB + D1_A - 2.0e0*D2_B)*(RAB + D1_A - 2.0e0*D2_B) + a_1_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB - D1_A - 2.0e0*D2_B)*(RAB - D1_A - 2.0e0*D2_B) + a_1_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + D1_A + 2.0e0*D2_B)*(RAB + D1_A + 2.0e0*D2_B) + a_1_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB - D1_A + 2.0e0*D2_B)*(RAB - D1_A + 2.0e0*D2_B) + a_1_2_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB + D1_A)*(RAB + D1_A) + a_1_2_quad))
		- 0.25e0*(1.0e0 / sqrt((RAB - D1_A)*(RAB - D1_A) + a_1_2_quad));

	return colchete;
}


double FourCenterIntegrals::miz_miz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D1_A = Params::get_double(mol.atom_name[A], "D1");
	double ro_1_A = Params::get_double(mol.atom_name[A], "ro_1");
	double D1_B = Params::get_double(mol.atom_name[B], "D1");
	double ro_1_B = Params::get_double(mol.atom_name[B], "ro_1");

	double a_1_1_quad = (ro_1_A + ro_1_B)*(ro_1_A + ro_1_B);

	double colchete = + 0.25e0*(1.0e0 / sqrt((RAB + D1_A - D1_B)*(RAB + D1_A - D1_B) + a_1_1_quad))
		- 0.25e0*(1.0e0 / sqrt((RAB + D1_A + D1_B)*(RAB + D1_A + D1_B) + a_1_1_quad))
		- 0.25e0*(1.0e0 / sqrt((RAB - D1_A - D1_B)*(RAB - D1_A - D1_B) + a_1_1_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB - D1_A + D1_B)*(RAB - D1_A + D1_B) + a_1_1_quad));

	return colchete;
}


double FourCenterIntegrals::mipi_Qpiz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D1_A = Params::get_double(mol.atom_name[A], "D1");
	double ro_1_A = Params::get_double(mol.atom_name[A], "ro_1");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_1_2_quad = (ro_1_A + ro_2_B)*(ro_1_A + ro_2_B);

	double colchete = -0.25e0*(1.0e0 / sqrt((RAB - D2_B)*(RAB - D2_B) + (D1_A - D2_B)*(D1_A - D2_B) + a_1_2_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB - D2_B)*(RAB - D2_B) + (D1_A + D2_B)*(D1_A + D2_B) + a_1_2_quad))
		+ 0.25e0*(1.0e0 / sqrt((RAB + D2_B)*(RAB + D2_B) + (D1_A - D2_B)*(D1_A - D2_B) + a_1_2_quad))
		- 0.25e0*(1.0e0 / sqrt((RAB + D2_B)*(RAB + D2_B) + (D1_A + D2_B)*(D1_A + D2_B) + a_1_2_quad));

	return colchete;
}


double FourCenterIntegrals::Qpiz_Qpiz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = Params::get_double(mol.atom_name[A], "D2");
	double ro_2_A = Params::get_double(mol.atom_name[A], "ro_2");
	double D2_B = Params::get_double(mol.atom_name[B], "D2");
	double ro_2_B = Params::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete = 0.125e0*(1.0e0 / sqrt((RAB + D2_A - D2_B)*(RAB + D2_A - D2_B) + (D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + D2_A - D2_B)*(RAB + D2_A - D2_B) + (D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB + D2_A + D2_B)*(RAB + D2_A + D2_B) + (D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB + D2_A + D2_B)*(RAB + D2_A + D2_B) + (D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB - D2_A - D2_B)*(RAB - D2_A - D2_B) + (D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB - D2_A - D2_B)*(RAB - D2_A - D2_B) + (D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))
		+ 0.125e0*(1.0e0 / sqrt((RAB - D2_A + D2_B)*(RAB - D2_A + D2_B) + (D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		- 0.125e0*(1.0e0 / sqrt((RAB - D2_A + D2_B)*(RAB - D2_A + D2_B) + (D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))		
		;
	return colchete;
}



