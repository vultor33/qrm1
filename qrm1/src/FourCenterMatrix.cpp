#include "FourCenterMatrix.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>

#include "FourCenterIntegrals.h"
#include "Molecule.h"
#include "MatrixRotation.h"
#include "Params.h"
#include "Coordstructs.h"

using namespace std;

void FourCenterMatrix::build_matrix_integrals(Molecule &mol)
{

	notANumber = false;
	int size = 1; // all spxpypz is defined even for H atoms.
	pMol_ = &mol;

	matrix_integrals.resize(mol.number_of_atoms);
	for(int i=0; i<mol.number_of_atoms;i++)
	{
		matrix_integrals[i].resize(mol.number_of_atoms);
		for(int j=0; j<mol.number_of_atoms; j++)
		{
			matrix_integrals[i][j].f_center.resize(size);
			for(int k=0; k<size; k++)
			{
				matrix_integrals[i][j].f_center[k].resize(size);
				for(int l=0; l<size; l++)
				{
					matrix_integrals[i][j].f_center[k][l].resize(size);
					for(int m=0; m<size; m++)
					{
						matrix_integrals[i][j].f_center[k][l][m].resize(size);
					}
				}
			}
		}
	}

	double value;
	bool not_done;

	//int mit,nit,lambdat,sigmat;

	for(int A=0; A<mol.number_of_atoms; A++)
	{
		int n_base_A = Params::get_int(mol.atom_name[A],"base_number")-1;

		for(int B=0; B<mol.number_of_atoms; B++)
		{
			int n_base_B = Params::get_int(mol.atom_name[B],"base_number")-1;

			not_done = true;
			for(int mi=0; mi<size; mi++)
			{
				if(mi>n_base_A)
				{
					value = 0.0e0;
					not_done = false;
				}
				else
				{
					not_done = true;
				}
				for(int ni=mi; ni<size; ni++)
				{
					if(ni>n_base_A)
					{
						value = 0.0e0;
						not_done = false;
					}
					else
					{
						not_done = true;
					}
					for(int lambda=0; lambda<size; lambda++)
					{	

						if(lambda>n_base_B)
						{
							value = 0.0e0;
							not_done = false;
						}
						else
						{
							not_done = true;
						}
						for(int sigma=lambda; sigma<size; sigma++)
						{					
							if(sigma>n_base_B)
							{
								value = 0.0e0;
								not_done = false;
							}
							else
							{
								not_done = true;
							}
							
							if((A==B)&&
								((lambda>n_base_A)||(sigma>n_base_A))&&
								((mi>n_base_B)||(ni>n_base_B)))
							{
								value = 0.0e0;
								not_done = false;
							}
							else
							{
								not_done = true;
							}

						/*	
							if(not_done)
							{
								value = integrals_.get_four_center_integral(mol,A,B,mi,ni,lambda,sigma);

							}
							//   cout << "mi:  " << mit << "  ni  " << nit << "  lambda  " << lambdat << "  sigma  " << sigmat
							//	   << "  value  " << value << endl;

							mit = transforma(mi);
							nit = transforma(ni);
							lambdat = transforma(lambda);
							sigmat = transforma(sigma);

								matrix_integrals[A][B].f_center[mit][nit][lambdat][sigmat] = value;
								matrix_integrals[A][B].f_center[mit][nit][sigmat][lambdat] = value;
								matrix_integrals[A][B].f_center[nit][mit][lambdat][sigmat] = value;
								matrix_integrals[A][B].f_center[nit][mit][sigmat][lambdat] = value;
						}
						*/

							if(not_done)
							{
								value = integrals_.get_four_center_integral(mol,A,B,mi,ni,lambda,sigma);
								if (isnan(value))
								{
									notANumber = true;
									value = 0.0e0;
								}
							}
							//   cout << "mi:  " << mi << "  ni  " << ni << "  lambda  " << lambda << "  sigma  " << sigma
							//	   << "  value  " << value << endl;

								matrix_integrals[A][B].f_center[mi][ni][lambda][sigma] = value;
								matrix_integrals[A][B].f_center[mi][ni][sigma][lambda] = value;
								matrix_integrals[A][B].f_center[ni][mi][lambda][sigma] = value;
								matrix_integrals[A][B].f_center[ni][mi][sigma][lambda] = value;
						}


					}
				}
			}

			//cin.get();

			if(A!=B)
			{
//				four_center_rotation(mol,A,B);
			}
		}
	}
}

int FourCenterMatrix::transforma(int base)
{
	switch(base)
	{
	case 0:
		return 0;
		break;
	case 1:
		return 3;
		break;
	case 2:
		return 1;
		break;
	case 3:
		return 2;
	default:
		return -1;
	}
}


void FourCenterMatrix::four_center_rotation(const Molecule &mol, int A, int B)
{
	int size = 4;
	vector< vector< vector< vector<double> > > > old_matrix = matrix_integrals[A][B].f_center;
	vector< vector<double> > rotation_matrix = m_rotation_.CalcRotatingMatrix(mol,A,B);
	//vector<vector<double>> rotation_matrix = m_rotation_.CalcRotatingMatrix_over(mol,A,B);

	fstream file;
	file.open("fred-debug-teste.txt", fstream::out);
	file << setiosflags(ios::fixed) << setprecision(6) << setw(12);

	for(int mi=0; mi<size; mi++){
		for(int ni=0; ni<size; ni++){
			for(int lambda=0; lambda<size; lambda++){
				for(int sigma=0; sigma<size; sigma++){
					matrix_integrals[A][B].f_center[mi][ni][lambda][sigma] = 0.0;
					for(int i=0; i<size; i++){
						for(int j=0; j<size; j++){
							for(int k=0; k<size; k++){
								for(int l=0; l<size; l++){
									matrix_integrals[A][B].f_center[mi][ni][lambda][sigma] +=
										old_matrix[i][j][k][l] 
										*rotation_matrix[mi][i] 
										*rotation_matrix[ni][j] 
										*rotation_matrix[lambda][k] 
										*rotation_matrix[sigma][l];
								}
							}
						}
					}
					file << " mi ni lambda sigma:   " << mi << "   " << ni << "   " << lambda << "   " << sigma << "  |   " <<  
						old_matrix[mi][ni][lambda][sigma] << "   r    " << matrix_integrals[A][B].f_center[mi][ni][lambda][sigma] << endl;
				}
			}
		}
	}
	file.close();
}


/*
double FourCenterMatrix::get_four_center(int A, int B, int mi, int ni, int lambda, int sigma)
{
	return matrix_integrals[A][B].f_center[mi][ni][lambda][sigma];
}
*/

double FourCenterMatrix::get_four_center(int A, int B, int mi, int ni, int lambda, int sigma, bool isExchange)
{
	double q = Params::get_double(pMol_->atom_name[A], "qalfa");
	if (
		(pMol_->scfMethod == "RHF") ||
		(pMol_->scfMethod == "UHF") ||
		(abs(q - 1.0e0) < 1.0e-6) ||
		(!isExchange)
		)
	{
		return matrix_integrals[A][B].f_center[mi][ni][lambda][sigma];
	}
	else
	{
		double R_distancia;
		if (A == B)
			R_distancia = 0.0e0;
		else
			R_distancia = pMol_->distancia_i_j[A][B];

		double ro_0_A = Params::get_double(pMol_->atom_name[A], "ro_0");
		double ro_0_B = Params::get_double(pMol_->atom_name[B], "ro_0");

		double integral_AB;
		double a00square = (ro_0_A + ro_0_B)*(ro_0_A + ro_0_B);
		double gss = Params::get_double(pMol_->atom_name[B], "gss");
		double linearAdjsutX = pow(1.5e0 - 0.5e0*pow(a00square, 1 - q), (1 / (1 - q)));
		
		double x = R_distancia*R_distancia + a00square;
		integral_AB = pow(1.5e0 - 0.5e0*pow(x, 1 - q), (1 / (1 - q)))
			-linearAdjsutX + gss;
		if (isnan(integral_AB))
		{
			notANumber = true;
			return 0.0e0;
		}
		else
			return integral_AB;
	}
}


