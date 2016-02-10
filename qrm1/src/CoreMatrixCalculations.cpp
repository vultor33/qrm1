#include "CoreMatrixCalculations.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "Params.h"

using namespace std;

CoreMatrixCalculations::CoreMatrixCalculations(Molecule &mol, FourCenterMatrix &four_center_in)
:four_center_(four_center_in)
{
	notANumber = false;
	build_fock_matrix(mol);
}

void CoreMatrixCalculations::build_fock_matrix(Molecule &mol)
{
	set_fock_matrix_dimension(mol);
	ressonance_AB_.start_overlaps(mol);

	int i_fock = -1;
	int j_fock = -1;
	int auxni;
	for (int A = 0; A < mol.number_of_atoms; A++)
	{
		int n_bases_A = Params::get_int(mol.atom_name[A], "base_number");
		for (int mi = 0; mi < n_bases_A; mi++)
		{
			i_fock++;
			mol.i_fock_base_line[A][mi] = i_fock;
			j_fock = i_fock;
			double fockElement;
			for (int B = A; B < mol.number_of_atoms; B++)
			{
				int n_bases_B = Params::get_int(mol.atom_name[B], "base_number");
				if (B == A)
				{
					auxni = mi;
				}
				else
				{
					auxni = 0;
				}
				for (int ni = auxni; ni < n_bases_B; ni++)
				{
					if ((A == B) && (mi == ni))
					{
						fock_matrix[i_fock][j_fock] = mi_mi_integral(mol, A, mi);
					}
					else if ((A == B))
					{
						fock_matrix[i_fock][j_fock] = mi_ni_integral(mol, A, mi,ni);
					}
					else
					{
						fockElement = ressonance_AB_.get_overlap_AB(A,B,mi,ni);// diatomic positions not fock positions
						if(isnan(fockElement))
						{
							notANumber = true;
							fockElement = 0.0e0;
						}
						fock_matrix[i_fock][j_fock] = fockElement;
					}
					//cout << "A: " << A << "B: " << B << "mi: " << mi << "ni: " << ni << "  i:  " << i_fock << "  j:  " << j_fock << endl;
					j_fock++;
				}
			}
		}
	}
}

void CoreMatrixCalculations::set_fock_matrix_dimension(Molecule &mol)
{
	int num_atomos = mol.number_of_atoms;
	mol.i_fock_base_line.resize(num_atomos);

	int dimensao_fock = 0;
	for (int i = 0; i < num_atomos; i++)
	{
		dimensao_fock += Params::get_int(mol.atom_name[i],"base_number");
		mol.i_fock_base_line[i].resize(Params::get_int(mol.atom_name[i], "base_number"));
	}
	fock_matrix.resize(dimensao_fock);
	for (int i = 0; i < dimensao_fock; i++)
	{
		fock_matrix[i].resize(dimensao_fock);
	}

}


double CoreMatrixCalculations::mi_mi_integral(const Molecule &mol, int A, int mi)
{
	double H_mi_mi;
	
	if (mi == 0)
	{
		H_mi_mi = Params::get_double(mol.atom_name[A], "uss");
	}
	else
	{
		H_mi_mi = Params::get_double(mol.atom_name[A], "upp");
	}

	for (int B = 0; B < mol.number_of_atoms; B++)
	{
		if (B != A)
		{
			H_mi_mi += mi_ni_potential(mol, A, B, mi,  mi);
		}
	}

	return H_mi_mi;
}

double CoreMatrixCalculations::mi_ni_integral(const Molecule &mol, int A, int mi, int ni)
{
	double H_mi_ni = 0.0e0;
	for (int B = 0; B < mol.number_of_atoms; B++)
	{
		if (B != A)
		{
			// so os orbitais s;
			H_mi_ni += mi_ni_potential(mol, A, B, mi, ni);
		}
	}

	return H_mi_ni;
}

double CoreMatrixCalculations::mi_ni_potential(const Molecule &mol, int A, int B, int mi, int ni)
{
	double potential = -((double)Params::get_int(mol.atom_name[B], "atomic_charge"))*
		four_center_.get_four_center(A,B,mi,ni,0,0);

	return potential;
}

void CoreMatrixCalculations::print_core_matrix()
{
	int n_bases = fock_matrix.size();
	cout << "Matriz de fock: " << endl;
	for (int i = 0; i < n_bases;i++)
	{
		for (int j = 0; j < n_bases;j++)
		{
			cout << fock_matrix[i][j] * Params::hartree_ev << "  ";
		}
		cout << endl;
	}
}

