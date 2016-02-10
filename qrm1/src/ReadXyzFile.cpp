#include "ReadXyzFile.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "Coordstructs.h"

using namespace std;

void ReadXyzFile::abrelexyz(string s)
{
	arqxyz.open(s.c_str());
	if(!arqxyz.is_open())
	{
		cout << "ERRO NA CLASSE: lexyz" << endl
			<< " O arquivo nao abriu corretamente" << endl;
		cin.get();
		exit(1);
	}

	arqxyz >> natoms;
	arqxyz.ignore(80,'\n');
	getline (arqxyz,flavor);
}

void ReadXyzFile::leiaCoordXYZ(vector<CoordXYZ> &molecula)
{
	for(int i=0;i<natoms;i++)
	{
		arqxyz >> molecula[i].atomlabel 
			>> molecula[i].x 
			>> molecula[i].y 
			>> molecula[i].z;
	}
	arqxyz.close();
}
