#include "ZmatoCart.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "Matrixop.h"
#include "Coordstructs.h"
#include "Params.h"

using namespace std;

void ZmatoCart::ztocart(vector<CoordZMAT> MolZ)
{

	int natoms = MolZ.size();
	if(natoms>0)
	{
		MolXYZ.resize(natoms);
	}
	else
	{
		cout << "Erro na rotina que converte coordenadas"
			<< endl
			<< "numero de atomos invalido"
			<< endl;
		cin.get();
		exit(2);
	}

	if(natoms==1)
	{
		ztocart1atomo(MolZ);
	}
	else if(natoms==2)
	{
		ztocart1atomo(MolZ);
		ztocart2atomos(MolZ);
	}
	else if(natoms==3)
	{
		ztocart1atomo(MolZ);
		ztocart2atomos(MolZ);
		ztocart3atomos(MolZ);
	}
	else if(natoms>3)
	{
		ztocart1atomo(MolZ);
		ztocart2atomos(MolZ);
		ztocart3atomos(MolZ);
		for(int i=3;i<natoms;i++)
		{
			ztocartiatomo(MolZ,i);
		}

	}
}

void ZmatoCart::ztocart1atomo(vector<CoordZMAT> MolZ)
{
	MolXYZ[0].atomlabel = MolZ[0].atomlabel;
	MolXYZ[0].x=0.0e0;
	MolXYZ[0].y=0.0e0;
	MolXYZ[0].z=0.0e0;
}

void ZmatoCart::ztocart2atomos(vector<CoordZMAT> MolZ)
{
	MolXYZ[1].atomlabel = MolZ[1].atomlabel;
	MolXYZ[1].x= MolZ[1].dis;
	MolXYZ[1].y=0.0e0;
	MolXYZ[1].z=0.0e0;
}

void ZmatoCart::ztocart3atomos(vector<CoordZMAT> MolZ)
{// os dois conect estao preenchidos, porque uso so um?
	double angrad = MolZ[2].ang*(Params::pi / 180);
	double auxy = MolZ[2].dis*sin(angrad);
	double auxxA = MolXYZ[MolZ[2].conect[0] - 1].x;
	double auxx;
	if (MolZ[2].conect[0] == 1)
		auxx = MolZ[2].dis*cos(angrad) + auxxA;
	else
		auxx = -MolZ[2].dis*cos(angrad) + auxxA;

	MolXYZ[2].atomlabel = MolZ[2].atomlabel;
	MolXYZ[2].x=auxx;
	MolXYZ[2].y=auxy;
	MolXYZ[2].z=0.0e0;
}

void ZmatoCart::ztocartiatomo(vector<CoordZMAT> MolZ, int i)
{

	double auxdis = -MolZ[i].dis;
	double auxteta = -(MolZ[i].ang)*(Params::pi / 180);
	double auxfi = MolZ[i].die*(Params::pi / 180);
	vector<double> vecjref(3);
	vecjref[0]=auxdis*cos(auxteta);
	vecjref[1]=auxdis*cos(auxfi)*sin(auxteta);
	vecjref[2]=auxdis*sin(auxfi)*sin(auxteta);

// AB = A e BC = B nesse momento
	vector<double> vetorAB(3);
	vetorAB[0] = MolXYZ[MolZ[i].conect[2]-1].x;
	vetorAB[1] = MolXYZ[MolZ[i].conect[2]-1].y;
	vetorAB[2] = MolXYZ[MolZ[i].conect[2]-1].z;

	vector<double> vetorBC(3);
	vetorBC[0] = MolXYZ[MolZ[i].conect[1]-1].x;
	vetorBC[1] = MolXYZ[MolZ[i].conect[1]-1].y;
	vetorBC[2] = MolXYZ[MolZ[i].conect[1]-1].z;

	vector<double> vecC(3);
	vecC[0] = MolXYZ[MolZ[i].conect[0]-1].x;
	vecC[1] = MolXYZ[MolZ[i].conect[0]-1].y;
	vecC[2] = MolXYZ[MolZ[i].conect[0]-1].z;

	vetorAB = vecAB(vetorAB[0],vetorAB[1],vetorAB[2],vetorBC[0],vetorBC[1],vetorBC[2]);

	vetorBC = vecAB(vetorBC[0],vetorBC[1],vetorBC[2],vecC[0],vecC[1],vecC[2]);

	vetorBC = normalizar(vetorBC[0],vetorBC[1],vetorBC[2]);

	vector<double> nchapeu(3);
	nchapeu = produtvetor(vetorAB[0],vetorAB[1],vetorAB[2],vetorBC[0],vetorBC[1],vetorBC[2]);
	nchapeu = normalizar(nchapeu[0],nchapeu[1],nchapeu[2]);

	vector<double> COL2(3);
	COL2 = produtvetor(nchapeu[0],nchapeu[1],nchapeu[2],vetorBC[0],vetorBC[1],vetorBC[2]);

	vector<double> vecDfim(3);
	vecDfim = matrizvec(vetorBC,COL2,nchapeu,vecjref);

	vecDfim = somavec(vecDfim[0],vecDfim[1],vecDfim[2],vecC[0],vecC[1],vecC[2]);

	MolXYZ[i].atomlabel = MolZ[i].atomlabel;
	MolXYZ[i].x=vecDfim[0];
	MolXYZ[i].y=vecDfim[1];
	MolXYZ[i].z=vecDfim[2];

}
