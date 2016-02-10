#include "DlibLinearSystem.h"

#include <iostream>
#include <dlib/matrix.h>
#include <vector>

using namespace std;
using namespace dlib;

DlibLinearSystem::~DlibLinearSystem(){}

bool DlibLinearSystem::calculateBCoefficients(std::vector< std::vector<double> > &bMatrix)
{
	long int size = bMatrix.size();
	std::vector<double> bVector = getbVector(size);

	matrix<double> M;
	matrix<double> y;

	M.set_size(size, size);
	y.set_size(size, 1);

	for (int i = 0; i < size; i++)
	{
		y(i) = bVector[i];
		for (int j = 0; j < size; j++)
		{
			M(i, j) = bMatrix[i][j];
		}
	}
	if (abs(det(M) < 1.0e-8))
	{
		return false;
	}
	else
	{
		matrix<double> x = inv(M)*y;

		std::vector<double> coefficientsOut(size - 1);
		for (int i = 1; i < size; i++)
		{
			coefficientsOut[i - 1] = x(i);
		}
		return true;
	}
}

std::vector<double> DlibLinearSystem::getCoefficientsOut()
{
	return coefficientsOut;
}

std::vector<double> DlibLinearSystem::getbVector(int size)
{
	std::vector<double> bVector(size);
	bVector[0] = -1.0e0;
	for (int i = 1; i < size; i++)
	{
		bVector[i] = 0.0e0;
	}
	return bVector;
}

/*
void DlibLinearSystem::teste()
{
	double a = 4;

	matrix<double, 3, 1> y;
	matrix<double> M(3, 3);

	M = 54.2, a, 12.1,
		1, 2, 3,
		5.9, 0.05, 1;

	M(0, 0) = a;

	y = 2.0,
		3.0,
		4.0;

	matrix<double> x = inv(M)*y;

	cout << "x: \n" << x << endl;

	cout << "M*x - y: \n" << M*x - y << endl;

}
*/