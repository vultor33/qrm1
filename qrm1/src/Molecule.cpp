#include "Molecule.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "ReadXyzFile.h"
#include "ReadMopFile.h"
#include "Coordstructs.h"
#include "Params.h"
#include "ZmatoCart.h"
#include "PrintAll.h"

using namespace std;

Molecule::Molecule(){}

//construir um q-input
Molecule::Molecule(string arquivoxyz, int buildType)
{
	if (buildType == 0)
	{
		get_extension_from_name(arquivoxyz);
		// trocar para SWITCH
		if (archive_extension == 0)
		{
			print_on_screen_XYZ(arquivoxyz);
		}
		else
		{
			build_molecule_by_ZMAT(arquivoxyz);
		}
	}
	// build 2
	// labels Hs
	// numero de atomos
	// coordenadas xyz


		atom_name.resize(number_of_atoms);
		distancia_i_j.resize(number_of_atoms);
		MolXYZ.resize(number_of_atoms);
		int aux_n_eletrons = 0;
		for (int i = 0; i < number_of_atoms; i++)
		{
			distancia_i_j[i].resize(number_of_atoms);
			if(archive_extension == 0)
				atom_name[i] = MolXYZ[i].atomlabel;
			else
				atom_name[i] = MolZMAT[i].atomlabel;
			aux_n_eletrons += Params::get_int(atom_name[i], "number_of_electrons");
		}
		number_of_electrons = aux_n_eletrons;
		if ((number_of_electrons % 2) != 0)
		{
			//number_of_electrons--;
			//charge = 1;
		}
		number_of_electrons--;
		charge = 1;

}

Molecule::Molecule(int charge_in, std::vector<std::string> labels)
{
	number_of_atoms = labels.size();
	atom_name.resize(number_of_atoms);
	distancia_i_j.resize(number_of_atoms);
	MolXYZ.resize(number_of_atoms);
	int aux_n_eletrons = 0;
	for (int i = 0; i < number_of_atoms; i++)
	{
		distancia_i_j[i].resize(number_of_atoms);
		atom_name[i] = labels[i];
		aux_n_eletrons += Params::get_int(atom_name[i], "number_of_electrons");
	}
	number_of_electrons = aux_n_eletrons - charge_in;
	charge = charge_in;
}

Molecule::~Molecule(){}

void Molecule::print_on_screen_XYZ(string arq)
{
	ReadXyzFile lecoord;
	lecoord.abrelexyz(arq);
	number_of_atoms=lecoord.natoms;
	MolXYZ.resize(number_of_atoms);
	lecoord.leiaCoordXYZ(MolXYZ);
}
void Molecule::build_molecule_by_ZMAT(string arq)
{
	ReadMopFile lecoord;
	lecoord.leiatudo(arq);
	number_of_atoms=lecoord.natoms;
	MolZMAT.resize(number_of_atoms);
	lecoord.pegacoordZMAT(MolZMAT);
}

void Molecule::print_molecule_on_screen()
{
	// trocar para SWITCH
	if(archive_extension==0)
	{
		print_on_screen();
	}
	else
	{
		print_on_screen_ZMAT();
	}
}

void Molecule::print_on_screen()
{
	cout << number_of_atoms << endl;
	cout << " texto inutil" << endl;

	for (int i=0;i<number_of_atoms;i++)
	{
		cout << MolXYZ[i].atomlabel << "  "
			<< MolXYZ[i].x << "   "
			<< MolXYZ[i].y << "   "
			<< MolXYZ[i].z << endl;
	}

}

void Molecule::print_on_screen_ZMAT()
{
	cout << number_of_atoms << endl;
	cout << " coordenadas internas no formato do mopac" << endl;

	for (int i=0;i<number_of_atoms;i++)
	{
		cout << MolZMAT[i].atomlabel 
			<< "   "
			<< MolZMAT[i].dis
			<< "   "
			<< MolZMAT[i].flagconec[0]
			<< "   "
			<< MolZMAT[i].ang
			<< "   "
			<< MolZMAT[i].flagconec[1]
			<< "   "
			<< MolZMAT[i].die
			<< "   "
			<< MolZMAT[i].flagconec[2]
			<< "   "
			<< MolZMAT[i].conect[0]
			<< "   "
			<< MolZMAT[i].conect[1]
			<< "   "
			<< MolZMAT[i].conect[2]
			<< endl;
	}

}

void Molecule::get_extension_from_name(string s)
{
	string extensaostring = s.substr(s.find_last_of(".") + 1);

//	char extensaostring[4];
//	int len = s.size();
//	int ultimo = s.copy(extensaostring,3,len-3);
//	extensaostring[ultimo]='\0';

	string xyz="xyz";

	if(extensaostring==xyz)
	{
		archive_extension=0;
	}
	else
	{
		archive_extension=1;

	}

}

void Molecule::build_atoms()
{
	ZmatoCart xyz;
	if (MolZMAT.size() != 0)
	{
		xyz.ztocart(MolZMAT);
		this->MolXYZ = xyz.MolXYZ;
	}
	for (int i = 0; i < number_of_atoms; i++)
	{
		for (int j = 0; j < number_of_atoms; j++)
		{
			if (i != j)
			{
				distancia_i_j[i][j] = Params::angs_bohr*
					sqrt((MolXYZ[i].x - MolXYZ[j].x)*(MolXYZ[i].x - MolXYZ[j].x)
					+ (MolXYZ[i].y - MolXYZ[j].y)*(MolXYZ[i].y - MolXYZ[j].y)
					+ (MolXYZ[i].z - MolXYZ[j].z)*(MolXYZ[i].z - MolXYZ[j].z));

			}
			else
			{
				distancia_i_j[i][j] = 0.0e0;
			}
		}
	}

	pPrintLog_->printMolecule(MolZMAT, MolXYZ, number_of_atoms);
}

void Molecule::buildXyzAtoms(vector<double> &coord)
{
	//coordinates: i   natoms+i   2*natoms+i
	for (int k = 0; k < number_of_atoms; k++)
	{
		MolXYZ[k].x = coord[k];
		MolXYZ[k].y = coord[k+number_of_atoms];
		MolXYZ[k].z = coord[k+2*number_of_atoms];
	}

	for (int i = 0; i < number_of_atoms; i++)
	{
		for (int j = 0; j < number_of_atoms; j++)
		{
			if (i != j)
			{
				distancia_i_j[i][j] = Params::angs_bohr*
					sqrt((MolXYZ[i].x - MolXYZ[j].x)*(MolXYZ[i].x - MolXYZ[j].x)
					+ (MolXYZ[i].y - MolXYZ[j].y)*(MolXYZ[i].y - MolXYZ[j].y)
					+ (MolXYZ[i].z - MolXYZ[j].z)*(MolXYZ[i].z - MolXYZ[j].z));
			}
			else
			{
				distancia_i_j[i][j] = 0.0e0;
			}
		}
	}

	pPrintLog_->printMolecule(MolZMAT, MolXYZ, number_of_atoms);

}


void Molecule::buildXyzAtoms(vector<CoordXYZ> &coord)
{
	MolXYZ = coord;
	for (int i = 0; i < number_of_atoms; i++)
	{
		for (int j = 0; j < number_of_atoms; j++)
		{
			if (i != j)
			{
				distancia_i_j[i][j] = Params::angs_bohr*
					sqrt((MolXYZ[i].x - MolXYZ[j].x)*(MolXYZ[i].x - MolXYZ[j].x)
					+ (MolXYZ[i].y - MolXYZ[j].y)*(MolXYZ[i].y - MolXYZ[j].y)
					+ (MolXYZ[i].z - MolXYZ[j].z)*(MolXYZ[i].z - MolXYZ[j].z));
			}
			else
			{
				distancia_i_j[i][j] = 0.0e0;
			}
		}
	}
	pPrintLog_->printMolecule(MolZMAT, MolXYZ, number_of_atoms);
}



/*
// apagar e substituir
void Molecule::build_atoms()
{
mol_atoms.resize(number_of_atoms);
int aux_n_eletrons = 0;
for (int i = 0; i < number_of_atoms; i++)
{
mol_atoms[i].build_atom(MolZMAT, i);
aux_n_eletrons += mol_atoms[i].num_eletrons;

}
number_of_electrons = aux_n_eletrons;

build_atoms(0);
}
*/

