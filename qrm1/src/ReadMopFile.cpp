#include "ReadMopFile.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include "Coordstructs.h"

using namespace std;

void ReadMopFile::leiatudo(string s)
{
	arqmop.open(s.c_str());
	if(!arqmop.is_open())
	{
		cout << "ERRO NA CLASSE: leiamop" << endl
			<< " O arquivo nao abriu corretamente" << endl;
		cin.get();
		exit(1);
	}

	string auxlinha1;
	getline (arqmop,auxlinha1);
	todaslinhas.push_back(auxlinha1);
	string auxlinha2;
	getline (arqmop,auxlinha2);
	todaslinhas.push_back(auxlinha2);
	string auxlinha3;
	getline (arqmop,auxlinha3);
	todaslinhas.push_back(auxlinha3);

	int naux;
	bool fim=1;
	int i=0;
	natoms=0;

	while(fim)
	{
		string auxlinha;
		getline(arqmop,auxlinha);
		naux = auxlinha.size();
		if(naux==0)
		{
			break;
		}
		natoms++;
		todaslinhas.push_back(auxlinha);
	}

	arqmop.close();
}


void ReadMopFile::pegacoordZMAT(vector<CoordZMAT> &molecula)
{
	for(int i=0;i<natoms;i++)
	{
		stringstream linha;
		linha << todaslinhas[i+3];

		linha >> molecula[i].atomlabel 
			>> molecula[i].dis
			>> molecula[i].flagconec[0]
			>> molecula[i].ang
			>> molecula[i].flagconec[1]
			>> molecula[i].die
			>> molecula[i].flagconec[2]
			>> molecula[i].conect[0]
			>> molecula[i].conect[1]
			>> molecula[i].conect[2];
	}
}




