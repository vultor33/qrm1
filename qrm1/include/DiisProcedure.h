#ifndef DIISPROCEDURE_H
#define DIISPROCEDURE_H

#include <vector>

#include "DlibLinearSystem.h"
#include "MatrixDiagonalization.h"
#include "PrintAll.h"

class DiisProcedure
{
public:
	DiisProcedure(
		PrintAll * pPrintLog_in_,
		std::vector< std::vector<double> > &fockMatrix_in,
		std::vector< std::vector<double> > &densityMatrix_in,
		int maxMatrix_in = 5
		);
	~DiisProcedure();

	void diisCalculation();
	bool diisOn;

private:
	int nMatrix;
	int maxMatrix;
	int scfIteration;
	void buildBMatrix();
	double bMatrixElement(int i, int j);
	void getNewFockMatrix();
	double maxError(std::vector< std::vector<double> > &error);
	void eraseOldMatrix();
	std::vector< std::vector<double> > &fockMatrix;
	std::vector< std::vector<double> > &densityMatrix;

	std::vector< std::vector<double> > matrix_product(const std::vector< std::vector<double> > &matrix_A, const std::vector< std::vector<double> > &matrix_B);
	std::vector< std::vector< std::vector<double> > > errorMatrix;
	std::vector< std::vector<double> > bMatrix;
	std::vector< std::vector< std::vector<double> > > oldFockMatrix;

	DlibLinearSystem solveB_;
	MatrixDiagonalization matrDiag_;
	PrintAll * pPrintLog_;

};


#endif