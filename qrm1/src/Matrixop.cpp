#include "Matrixop.h"
#include <iostream>
#include <cmath>
#include <vector>

#include "Coordstructs.h"

using namespace std;

double Matrixop::produtescalar(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return x1*x2+y1*y2+z1*z2;
}

vector<double> Matrixop::produtvetor(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> prodvec(3);
	prodvec[0] = y1*z2-y2*z1;
	prodvec[1] = z1*x2-z2*x1;
	prodvec[2] = x1*y2-x2*y1;

	return prodvec;
}

vector<double> Matrixop::normalizar(double x1, double y1, double z1)
{
	double norma2 = produtescalar(x1, y1, z1, x1, y1, z1);
	double norma = sqrt(norma2);
	vector<double> vecnormal(3);
	vecnormal[0] = x1/norma;
	vecnormal[1] = y1/norma;
	vecnormal[2] = z1/norma;

	return vecnormal;
}

vector<double> Matrixop::somavec(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> somavec(3);
	somavec[0] = x1+x2;
	somavec[1] = y1+y2;
	somavec[2] = z1+z2;

	return somavec;
}

vector<double> Matrixop::vecAB(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> vecAB(3);
	vecAB[0] = -x1+x2;
	vecAB[1] = -y1+y2;
	vecAB[2] = -z1+z2;

	return vecAB;
}

// os vetores sao as colunas da matriz
vector<double> Matrixop::matrizvec(vector<double> COL1, vector<double> COL2, vector<double> COL3, vector<double> VECX)
{
	vector<double> VECfim(3);

	VECfim[0] = COL1[0]*VECX[0]+COL2[0]*VECX[1]+COL3[0]*VECX[2];
	VECfim[1] = COL1[1]*VECX[0]+COL2[1]*VECX[1]+COL3[1]*VECX[2];
	VECfim[2] = COL1[2]*VECX[0]+COL2[2]*VECX[1]+COL3[2]*VECX[2];

	return VECfim;
}

vector< vector<double> > Matrixop::matrix_matrix(vector< vector<double> > matA, vector< vector<double> > matB)
{
	vector< vector<double> > matC;
	int size = matA.size();
	matC.resize(size);
	for(int i=0; i<size;i++)
	{
		matC[i].resize(size);
		for(int j=0; j<size;j++)
		{
			matC[i][j]=0.0e0;
			for(int k=0;k<size;k++)
			{
				matC[i][j] +=matA[i][k]*matB[k][j];
			}
		}
	}

	return matC;
}

