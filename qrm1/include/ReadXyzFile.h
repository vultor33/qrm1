#ifndef READXYZFILE_H
#define READXYZFILE_H

#include <string>
#include <fstream>
#include <vector>

#include "Coordstructs.h"

class ReadXyzFile
{
	std::string flavor;
	std::ifstream arqxyz;
public:
	int natoms;
	void abrelexyz(std::string);
	void leiaCoordXYZ(std::vector<CoordXYZ> &);
};

#endif


