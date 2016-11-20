#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <string>
#include <vector>

#include "ReadXyzFile.h"
#include "ReadMopFile.h"
#include "Coordstructs.h"
#include "PrintAll.h"

class Molecule
{
public:
	Molecule();
	Molecule(std::string, int buildType=0);
	Molecule(int charge_in, std::vector<std::string> labels);
	~Molecule();
	int number_of_atoms;
	int number_of_electrons;
	std::string scfMethod;
	std::vector<CoordXYZ> MolXYZ;
	std::vector<CoordZMAT> MolZMAT;

	void print_molecule_on_screen();
	inline void setPrintLog(PrintAll & printLog_in_){ pPrintLog_ = &printLog_in_; }

	//funcoes do novo - privadas;
	std::vector< std::vector<double> > distancia_i_j;
	std::vector<std::string> atom_name;
	std::vector< std::vector<int> > i_fock_base_line;

	// publicas
	void build_atoms();
	void buildXyzAtoms(std::vector<double> &coord);
	void buildXyzAtoms(std::vector<CoordXYZ> &coord);
	void print_on_screen();

	int getCharge() { return charge; }

private:
	int charge;
	int archive_extension;
	void print_on_screen_ZMAT();
	void print_on_screen_XYZ(std::string);
	void build_molecule_by_ZMAT(std::string);
	void get_extension_from_name(std::string);
	
	PrintAll * pPrintLog_;

};

#endif


// extensao
// 0 -> .xyz
// 1 -> .mop
