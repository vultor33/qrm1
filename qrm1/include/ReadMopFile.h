// ATENCAO - PROBLEMA NESSA ROTINA
// ELA NAO ESTA LENDO O NUMERO DE ATOMOS COMO DEVERIA
// E SIM CRIANDO VETORES MUITO GRANDES.

#ifndef READMOPFILE_H
#define READMOPFILE_H

#include <string>
#include <fstream>
#include <vector>

#include "Coordstructs.h"

class ReadMopFile
{
	std::string flavor[3];
	std::ifstream arqmop;
	std::vector<std::string> todaslinhas;
public:
	int natoms;
	void leiatudo(std::string);
	void pegacoordZMAT(std::vector<CoordZMAT> &);
};

#endif




