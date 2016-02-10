#ifndef POPLEOVELAPS_H
#define POPLEOVELAPS_H

#include <vector>

#include "Molecule.h"

class PopleOverlaps
{
private:
	double GetReducedOverlapAOs(int na, int nb, double alpha, double beta);
	double GetReducedOverlapAOs(int na, int la, int m, int nb, int lb, double alpha, double beta);
	double GetAuxiliaryA(int k, double rho);
	double GetAuxiliaryB(int k, double rho);
	double GetAuxiliaryD(int la, int lb, int m);
	int factorial(int n);
	int getM(int orbital);
	int getL(int orbital);
	void CalcDiatomicOverlapAOsInDiatomicFrame(const Molecule &mol,	int A, int B,  std::vector< std::vector<double> > &diatomicOverlapAOs);
	static const double Z[6][6][11];
	static const double Y[4][4][3][3][3][7][7];

public:
	std::vector< std::vector<double> > overlap_control(const Molecule &mol, int A, int B);
};



#endif