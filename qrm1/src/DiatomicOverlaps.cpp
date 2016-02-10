#include "DiatomicOverlaps.h"
#include<iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "Molecule.h"
#include "Params.h"
#include "CheckRessonanceType.h"
#include "Matrixop.h"
#include "MatrixRotation.h"
#include "PopleOverlaps.h"

using namespace std;


void DiatomicOverlaps::start_overlaps(const Molecule &mol)
{
	overlaps_AB.resize(mol.number_of_atoms);
	for(int i=0; i<mol.number_of_atoms; i++)
	{
		overlaps_AB[i].resize(mol.number_of_atoms);
	}

	for(int A=0; A<mol.number_of_atoms; A++)
	{
		for(int B=0; B<mol.number_of_atoms; B++)
		{
			if(A!=B)
			{

				overlaps_AB[A][B].bases = calc_overlap_AB(mol, A, B);
//				overlaps_AB[A][B].bases = popOver_.overlap_control(mol,A,B);
//				overlapToRessonance(mol,overlaps_AB[A][B].bases,A,B);

			}
		}
	}

}


void DiatomicOverlaps::overlapToRessonance(const Molecule &mol, vector< vector<double> > &over_AB, int A, int B)
{
	int n_basesA = Params::get_int(mol.atom_name[A],"base_number");
	int n_basesB = Params::get_int(mol.atom_name[B],"base_number");

//	utilidades util_;
//	cout << "overlap antes de rodar " << endl;
//	util_.print_matrix_on_file(over_AB);
//	cin.get();
//	util_.print_matrix_on_screen(over_AB);


	pople_rotation(mol, over_AB, A, B);


//	vector<vector<double>> temp_overlap = matrix_operation_.matrix_matrix(rotation_matrix,over_AB);
//	over_AB = matrix_operation_.matrix_matrix(temp_overlap,rotation_matrix);
//	cout << "overlap depois " << endl;	
//	util_.print_matrix_on_file(over_AB);
//	cout << "rotating " << endl;
//	util_.print_matrix_on_screen(rotation_matrix);
//	cin.get();

	int size = 4;
	for(int mi=0; mi<size;mi++)
	{
		for(int ni=0; ni<size;ni++)
		{
			if((n_basesA<(mi+1))||(n_basesB<(ni+1)))// always 4x4. If atom=H pxpypz is zero
			{
				over_AB[mi][ni]=0.0e0;
			}
			else
			{
				over_AB[mi][ni]*=0.5e0*(getBetaAtom(mol,mi,A)+getBetaAtom(mol,ni,B));
			}
		}
	}
}


void DiatomicOverlaps::pople_rotation(const Molecule &mol, vector< vector<double> > &over_AB, int A, int B)
{
	//vector<vector<double>> tmat = m_rotation_.CalcRotatingMatrix_over(mol,A,B);
	vector< vector<double> > tmat = m_rotation_.CalcRotatingMatrix(mol,A,B);
	vector< vector<double> > overTemp = over_AB;

	int size = tmat.size();
	for(int mi=0; mi<size; mi++)
	{
		for(int ni=0; ni<size; ni++)
		{
			overTemp[mi][ni] = 0.0e0;
			for(int a=0; a<size; a++)
			{
				for(int b=0; b<size; b++)
				{
					overTemp[mi][ni] += tmat[mi][a]*over_AB[a][b]*tmat[ni][b];
				}
			}
		}
	}

//	utilidades util_;
//	cout << "roda roda " << endl;	
//	util_.print_matrix_on_screen(overTemp);
/*
	over_AB[0][0] = overTemp[0][0]; //ss
	over_AB[0][1] = overTemp[0][3]; //spx
	over_AB[0][2] = overTemp[0][1]; //spy
	over_AB[0][3] = overTemp[0][2]; //spz
	over_AB[1][0] = overTemp[3][0]; //pxs
	over_AB[1][1] = overTemp[3][3]; //pxpx
	over_AB[1][2] = overTemp[3][1]; //pxpy
	over_AB[1][3] = overTemp[3][2]; //pxpz
	over_AB[2][0] = overTemp[1][0]; //pys
	over_AB[2][1] = overTemp[1][3]; //pypx
	over_AB[2][2] = overTemp[1][1]; //pypy
	over_AB[2][3] = overTemp[1][2]; //pypz
	over_AB[3][0] = overTemp[2][0]; //pzs
	over_AB[3][1] = overTemp[2][3]; //pzpx
	over_AB[3][2] = overTemp[2][1]; //pzpy
	over_AB[3][3] = overTemp[2][2]; //pzpz
*/
	over_AB = overTemp;
//	cout << "over ab novo " << endl;	
//	util_.print_matrix_on_screen(over_AB);
//	cin.get();
}






double DiatomicOverlaps::getBetaAtom(const Molecule &mol, int orbital, int Atom)
{
	switch(orbital)
	{
	case 0:
		return Params::get_double(mol.atom_name[Atom],"betas");
	case 1:
	case 2:
	case 3:
		return Params::get_double(mol.atom_name[Atom],"betap");
	default:
		cout << "erro em PopleOverlaps::getBetaAtom " << endl;
		exit(7);
	}

}


double DiatomicOverlaps::get_overlap_AB(int A, int B, int mi, int ni)
{
	return overlaps_AB[A][B].bases[mi][ni];
}








































vector< vector<double> > DiatomicOverlaps::calc_overlap_AB(const Molecule &mol, int A, int B)
{
	vector< vector<double> > over_AB;
	int n_basesA = Params::get_int(mol.atom_name[A],"base_number");
	int n_basesB = Params::get_int(mol.atom_name[B],"base_number");

	int size = 4;
	over_AB.resize(size);
	for(int mi=0; mi<size;mi++)
	{
		over_AB[mi].resize(size);
		for(int ni=0; ni<size;ni++)
		{
			if((n_basesA<(mi+1))||(n_basesB<(ni+1)))// always 4x4. If atom=H pxpypz is zero
			{
				over_AB[mi][ni]=0.0e0;
			}
			else
			{
				over_AB[mi][ni]=calc_overlap_AB_mini(mol,A,B,mi,ni);
			}
		}
	}
//	utilidades util_;
//	util_.print_matrix_on_screen(over_AB);
//	cin.get();

	vector< vector<double> > rotation_matrix = m_rotation_.CalcRotatingMatrix(mol,A,B);
	vector< vector<double> > temp_overlap = matrix_operation_.matrix_matrix(rotation_matrix,over_AB);
	over_AB = matrix_operation_.matrix_matrix(temp_overlap,rotation_matrix);

	return over_AB;

}

// calculando no metodo antigo
double DiatomicOverlaps::calc_overlap_AB_mini(const Molecule &mol, int A, int B, int mi, int ni)
{
	// A!=B for definition
	// todos os expoentes sao iguais, mas em um processo de parametrizacao isso nao seria
	// necessariamente verdadeiro. Isso pode dar problema em varias integrais.
	int q_number_A = Params::get_int(mol.atom_name[A], "q_number");
	int q_number_B = Params::get_int(mol.atom_name[B], "q_number");
	switch (check_.choose_ressonance_expression(q_number_A, q_number_B, mi, ni))
	{
	case 0:
		return 0.0e0;
	case 1:
		return overlap_1s_1s(mol, A, B);
	case 2:
		return overlap_1s_2s(mol, A, B);
	case 3:
		return overlap_1s_2s(mol, B, A);
	case 4:
		return overlap_1s_2pz(mol, A, B);
	case 5:
		return overlap_1s_2pz(mol, B, A);
	case 6:
		return overlap_2s_2s(mol, A, B);
	case 7:
		return overlap_2s_2pz(mol, A, B);
	case 8:
		return overlap_2s_2pz(mol, B, A);
	case 9:
		return overlap_2pz_2pz(mol, A, B);
	case 10:
		return overlap_2ppi_2ppi(mol, A, B);
	default:
		return 0.0e0;
		cout << "nao sei calcular essa integral - erro em ressonance_integral" << endl;
		exit(4);
	}
}

//A Study of Two-Center Integrals Useful in Calculations on Molecular Structure
// C. c. J. ROOTHAAN, 1951
double DiatomicOverlaps::overlap_1s_1s(const Molecule &mol, int A, int B)
{
	// so tem um tipo de integral 1s1s, e a entre hidrogenios
	double ro_expoente = Params::get_double(mol.atom_name[A], "expoents") * mol.distancia_i_j[A][B];

	double q = Params::get_double(mol.atom_name[A], "qover");
	double overlap;
	if (
		(Params::method[0] == 'q')&&
		(abs(1.0e0-q)>0.000001)
		)
	{
		if(ro_expoente > (1.0e0/(1.0e0 - q)) )
		{
			overlap = 0.0e0;
		}
		else
		{
			overlap = pow(1.0e0 - (1.0e0 - q)*ro_expoente, 1.0e0 / (1.0e0 - q))
			*(1.0e0 + ro_expoente + ro_expoente*ro_expoente / (3.0e0));
		}
	}
	else
	{
		overlap = (1.0e0 + ro_expoente + ro_expoente*ro_expoente / (3.0e0))*exp(-ro_expoente);
	}

	return overlap*Params::get_double(mol.atom_name[A], "betas");
}

//para parametrizacao o tau = 0 vai certamente acontecer
double DiatomicOverlaps::overlap_1s_2s(const Molecule &mol, int A, int B)
{
	double tau = (Params::get_double(mol.atom_name[A], "expoents") - Params::get_double(mol.atom_name[B], "expoents")) /
		(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoents"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoents"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = ( 1.0e0 - tau)*ro;

	if (ro == 0)
	{
		cout << "erro ro=0 consertar integrais " << endl;
		exit(5);
	}
	if (tau == 0)
	{
		cout << "erro tau=0 consertar integrais " << endl;
		exit(5);
	}
	double auxoverlap1 = sqrt(1.0e0 - tau*tau) / (sqrt(3.0e0)*tau*ro);
	double auxoverlap2 = -(1.0e0 - kapa)*
		(2.0e0*(1.0e0 + kapa)*(2.0e0 - 3.0e0*kapa) +
		(1.0e0 - 2.0e0*kapa)*ro_a)*exp(-ro_a);
	double auxoverlap3 = (1.0e0 + kapa)*
		(2.0e0*(1.0e0 - kapa)*(2.0e0 - 3.0e0*kapa) +
		4.0e0 * (1.0e0 - kapa)*ro_b + ro_b*ro_b)*exp(-ro_b);

	double overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);

	double H_mi_ni = overlap*0.5e0*(Params::get_double(mol.atom_name[A],"betas") +
		Params::get_double(mol.atom_name[B], "betas"));
	return H_mi_ni;
}

double DiatomicOverlaps::overlap_1s_2pz(const Molecule &mol, int A, int B)
{

	double tau = (Params::get_double(mol.atom_name[A], "expoents") - Params::get_double(mol.atom_name[B], "expoentp")) /
		(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoentp"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoentp"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = (1.0e0 - tau)*ro;

	if ((tau == 0) || (ro == 0))
	{
		cout << "erro tau=0 ou ro=0 consertar integrais " << endl;
		exit(5);
	}
	//atencao: essa divisao esta estranha
	double auxoverlap1 = sqrt((1.0e0 + tau) / (1.0e0 - tau))*
		(( 1.0e0) / (tau*ro*ro));
	double auxoverlap2 = -(1.0e0 - kapa)*(1.0e0 - kapa)*(
		6.0e0*(1.0e0 + kapa) * (1.0e0 + ro_a)
		+ 2.0e0 * ro_a*ro_a
		)*exp(-ro_a);
	double auxoverlap3 = (1.0e0 + kapa)*(
		6.0e0*(1.0e0 - kapa)*(1.0e0 - kapa)*(1.0e0 + ro_b) +
		4.0e0*(1.0e0 - kapa)*ro_b*ro_b + ro_b*ro_b*ro_b
		)*exp(-ro_b);

	double overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);

	double H_mi_ni = overlap*0.5e0*(
		Params::get_double(mol.atom_name[A], "betas") +
		Params::get_double(mol.atom_name[B], "betap"));

	return H_mi_ni;
}

double DiatomicOverlaps::overlap_2s_2s(const Molecule &mol, int A, int B)
{
	double tau = (Params::get_double(mol.atom_name[A], "expoents") - Params::get_double(mol.atom_name[B], "expoents")) /
		(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoents"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoents"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = (1.0e0 - tau)*ro;

	double overlap;

	if (ro == 0)
	{
		cout << "erro ro=0 consertar integrais " << endl;
		exit(5);
	}
	else if (tau == 0)
	{
		overlap = (1.0e0 + ro + (4.0e0 / 9.0e0)*ro*ro + (1.0e0 / 9.0e0)*ro*ro*ro
			+(1.0e0/45.0e0)*ro*ro*ro*ro)*exp(-ro);
	}
	else
	{
		double auxoverlap1 = sqrt(1.0e0 - tau*tau) / (3.0e0*tau*ro);
		double auxoverlap2 = -(1.0e0 - kapa)*(
			2.0e0*(1.0e0 + kapa)*(7.0e0 - 12.0e0*kapa*kapa) + 4.0e0*(1.0e0 + kapa)*(2.0e0 - 3.0e0*kapa)
			*ro_a + (1.0e0 - 2.0e0*kapa)*ro_a*ro_a)*exp(-ro_a);
		double auxoverlap3 = (1.0e0 + kapa)*(
			2.0e0*(1.0e0 - kapa)*(7.0e0 - 12.0e0*kapa*kapa) + 4.0e0*(1.0e0 - kapa)*(2.0e0 + 3.0e0*kapa)
			*ro_b + (1.0e0 + 2.0e0*kapa)*ro_b*ro_b)*exp(-ro_b);
		overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);
	}

	double H_mi_ni = overlap*0.5e0*(
		Params::get_double(mol.atom_name[A], "betas") +
		Params::get_double(mol.atom_name[B], "betas"));

	return H_mi_ni;

}

double DiatomicOverlaps::overlap_2s_2pz(const Molecule &mol, int A, int B)
{

	double tau = (Params::get_double(mol.atom_name[A], "expoents") - Params::get_double(mol.atom_name[B], "expoentp")) /
		(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoentp"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoents") + Params::get_double(mol.atom_name[B], "expoentp"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = (1.0e0 - tau)*ro;

	double overlap;

	if (ro == 0)
	{
		cout << "ro nao deveria ser zero - erro em ressonance 2s 2px" << endl;
		exit(5);
	}
	if (tau == 0)
	{
		double auxoverlap1 = (1.0e0 / (2.0e0*sqrt(3.0e0)))*ro;
		double auxoverlap2 = 1.0e0 + ro + (7.0e0 / 15.0e0)*ro*ro +
			(2.0e0 / 15.0e0)*ro*ro*ro;
		overlap = auxoverlap1*auxoverlap2*exp(-ro);
	}
	else
	{

		double auxoverlap1 = sqrt((1.0e0 + tau) / (1.0e0 - tau))*
			((1.0e0) / (sqrt(3.0e0)*tau*ro*ro));
		double auxoverlap2 = -(1.0e0 - kapa)*(1.0e0 - kapa)*(
			6.0e0*(1.0e0 + kapa)*(3.0e0 + 4.0e0*kapa)
			* (1.0e0 + ro_a) + 2.0e0 *(5.0e0 + 6.0e0*kapa)*ro_a*ro_a
			+ 2.0e0*ro_a*ro_a*ro_a)*exp(-ro_a);
		double auxoverlap3 = (1.0e0 + kapa)*(
			6.0e0*(1.0e0 - kapa)*(1.0e0 - kapa)*(3.0e0 + 4.0e0*kapa)*
			(1.0e0 + ro_b) +
			4.0e0*(1.0e0 - kapa)*(2.0e0 + 3.0e0*kapa)*ro_b*ro_b +
			(1.0e0 + 2.0e0*kapa)*ro_b*ro_b*ro_b
			)*exp(-ro_b);

		overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);
	}

	double H_mi_ni = overlap*0.5e0*(
		Params::get_double(mol.atom_name[A], "betas") +
		Params::get_double(mol.atom_name[B], "betap"));
	return H_mi_ni;
}


double DiatomicOverlaps::overlap_2pz_2pz(const Molecule &mol, int A, int B)
{
	double tau = (Params::get_double(mol.atom_name[A], "expoentp") - Params::get_double(mol.atom_name[B], "expoentp")) /
		(Params::get_double(mol.atom_name[A], "expoentp") + Params::get_double(mol.atom_name[B], "expoentp"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoentp") + Params::get_double(mol.atom_name[B], "expoentp"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = (1.0e0 - tau)*ro;

	double overlap;

	if (ro == 0)
	{
		cout << "ro nao deveria ser zero - erro em ressonance 2s 2px" << endl;
		exit(5);
	}
	if (tau == 0)
	{
		overlap = (-1.0e0 -ro -(1.0e0/5.0e0)*ro*ro + 
			(2.0e0/15.0e0)*ro*ro*ro + (1.0e0/15.0e0)*ro*ro*ro*ro)*exp(-ro);
	}
	else
	{

		double auxoverlap1 = 1.0e0 / (sqrt(1.0e0 - tau*tau)*tau*ro*ro*ro);
		double auxoverlap2 = -(1.0e0-kapa)*(1.0e0-kapa)*(
			24.0e0*(1.0e0 + kapa)*(1.0e0 + kapa)*(1.0e0+ro_a) +
			12.0e0*(1.0e0 + kapa)*ro_a*ro_a + 2.0e0*ro_a*ro_a*ro_a
			)*exp(-ro_a);
		double auxoverlap3 = (1.0e0+kapa)*(1.0e0+kapa)*(
			24.0e0*(1.0e0 - kapa)*(1.0e0 - kapa)*(1+ro_b)+
			12.0e0*(1.0e0 - kapa)*ro_b*ro_b + 2.0e0*ro_b*ro_b*ro_b
			)*exp(-ro_b);

		overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);
	}

	double H_mi_ni = overlap*0.5e0*(
		Params::get_double(mol.atom_name[A], "betap") +
		Params::get_double(mol.atom_name[B], "betap"));
	return H_mi_ni;

}


double DiatomicOverlaps::overlap_2ppi_2ppi(const Molecule &mol, int A, int B)
{
	double tau = (Params::get_double(mol.atom_name[A], "expoentp") - Params::get_double(mol.atom_name[B], "expoentp")) /
		(Params::get_double(mol.atom_name[A], "expoentp") + Params::get_double(mol.atom_name[B], "expoentp"));
	double zeta = 0.5e0*(Params::get_double(mol.atom_name[A], "expoentp") + Params::get_double(mol.atom_name[B], "expoentp"));
	double ro = zeta*mol.distancia_i_j[A][B];
	double kapa = 0.5e0*(tau + ((1.0e0) / tau));
	double ro_a = (1.0e0 + tau)*ro;
	double ro_b = (1.0e0 - tau)*ro;

	double overlap;

	if (ro == 0)
	{
		cout << "ro nao deveria ser zero - erro em ressonance 2s 2px" << endl;
		exit(5);
	}
	if (tau == 0)
	{
		overlap = (1.0e0 + ro + (2.0e0 / 5.0e0)*ro*ro + (1.0e0 / 15.0e0)*ro*ro*ro)
			*exp(-ro);
	}
	else
	{
		double auxoverlap1 = 1.0e0/(sqrt(1.0e0-tau*tau)*tau*ro*ro*ro);
		double auxoverlap2 = -(1.0e0 - kapa)*(1.0e0 - kapa)*(
			24.0e0*(1.0e0 + kapa)*(1.0e0 + kapa)*(1.0e0 + ro_a) +
			12.0e0*(1.0e0+kapa)*ro_a*ro_a + 2.0e0*ro_a*ro_a*ro_a
			)*exp(-ro_a);
		double auxoverlap3 = (1.0e0 + kapa)*(1.0e0 + kapa)*(
			24.0e0*(1.0e0 - kapa)*(1.0e0 - kapa)*(1.0e0+ro_b)+
			12.0e0*(1.0e0-kapa)*ro_b*ro_b + 2.0e0*ro_b*ro_b*ro_b
			)*exp(-ro_b);

		overlap = auxoverlap1*(auxoverlap2 + auxoverlap3);
	}

	double H_mi_ni = overlap*0.5e0*(
		Params::get_double(mol.atom_name[A], "betap") +
		Params::get_double(mol.atom_name[B], "betap"));
	return H_mi_ni;


}

