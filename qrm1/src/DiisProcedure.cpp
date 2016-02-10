#include "DiisProcedure.h"

#include <vector>

#include "DlibLinearSystem.h"
#include "PrintAll.h"

using namespace std;

DiisProcedure::DiisProcedure(
	PrintAll * pPrintLog_in_,
	vector< vector<double> > &fockMatrix_in,
	vector< vector<double> > &densityMatrix_in,
	int maxMatrix_in
	)
	:fockMatrix(fockMatrix_in), densityMatrix(densityMatrix_in)
{
	pPrintLog_ = pPrintLog_in_;
	diisOn = false;
	maxMatrix = maxMatrix_in;
	nMatrix = 0;
	scfIteration = 0;

}

DiisProcedure::~DiisProcedure(){}


//http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming
void DiisProcedure::diisCalculation()
{
	scfIteration++;
	if (scfIteration == 5)
	{
		pPrintLog_->printStartDiis();
		nMatrix = 0;
		diisOn = true;
	}

	if (diisOn)
	{
		matrDiag_.symmetrizeEntryMatrix(fockMatrix);
		vector< vector<double> > error = matrix_product(fockMatrix, densityMatrix);
		vector< vector<double> > errorTemp2 = matrix_product(densityMatrix, fockMatrix);
		int size = error.size();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				error[i][j] -= errorTemp2[i][j];
			}
		}
		double mError = maxError(error);
		pPrintLog_->printErrorDiis(mError);

		errorMatrix.push_back(error);
		oldFockMatrix.push_back(fockMatrix);
		buildBMatrix();
		eraseOldMatrix();
		pPrintLog_->printScfMatrix("diis error", error);
		getNewFockMatrix();
		pPrintLog_->printScfHeader(1);
	}
}

double DiisProcedure::maxError(vector< vector<double> > &error)
{
	double mError = error[0][0];
	int size = error.size();
	for(int i=0; i<size; i++)
	{
		for(int j=0;j<size; j++)
		{
			if(mError<error[i][j])
				mError = error[i][j];
		}
	}
	return mError;
}


void DiisProcedure::eraseOldMatrix()
{
	nMatrix++;
	if (nMatrix > 5)
	{
		errorMatrix.erase(errorMatrix.begin());
		oldFockMatrix.erase(oldFockMatrix.begin());
		bMatrix.erase(bMatrix.begin() + 1);
		for (int i = 0; i < (int)bMatrix[0].size(); i++)
		{
			bMatrix[i].erase(bMatrix[i].begin() + 1);
		}
	}
}

void DiisProcedure::getNewFockMatrix()
{
	if (solveB_.calculateBCoefficients(bMatrix))
	{
		vector<double> coefficients = solveB_.getCoefficientsOut();
		pPrintLog_->printScfMatrix("diis b coefficients", coefficients);
		int nMatr = errorMatrix.size();
		int size = errorMatrix[0].size();
		pPrintLog_->printScfMatrix("diis old fock", fockMatrix);
		for (int i = 0; i < size; i++)
		{
			for (int j = i; j < size; j++)
			{

				for (int A = 0; A < nMatr; A++)
				{
					fockMatrix[i][j] += coefficients[A] *
						oldFockMatrix[A][i][j];
				}
			}
		}
		matrDiag_.symmetrizeEntryMatrix(fockMatrix);
		pPrintLog_->printScfMatrix("diis new fock", fockMatrix);
	}
}


void DiisProcedure::buildBMatrix()
{
	// upper triangular
	if (bMatrix.size() == 0)
	{
		int bSize = errorMatrix.size() + 1;
		bMatrix.resize(bSize);
		for (int i = 0; i < bSize; i++)
		{
			bMatrix[i].resize(bSize);
		}
		bMatrix[0][0] = 0.0e0;
		bMatrix[0][1] = -1.0e0;
		bMatrix[1][1] = bMatrixElement(0, 0);
	}
	else
	{
		int bSize = errorMatrix.size() + 1;
		vector< vector<double> > bMatrixTemp;
		bMatrixTemp.resize(bSize);
		for (int i = 0; i < bSize; i++)
		{
			bMatrixTemp[i].resize(bSize);
		}

		for (int i = 0; i < bSize; i++)
		{
			for (int j = i; j < bSize; j++)
			{
				if ((i == 0) && (j == 0))
					bMatrixTemp[i][j] = 0.0e0;
				else if (i == 0)
					bMatrixTemp[i][j] = -1.0e0;
				else if (j == 0)
					bMatrixTemp[i][j] = -1.0e0;
				else if (j == (bSize - 1))
				{
					bMatrixTemp[i][j] = bMatrixElement(i - 1, j - 1);
				}
				else
					bMatrixTemp[i][j] = bMatrix[i][j];
			}
		}
		bMatrix = bMatrixTemp;
	}
	matrDiag_.symmetrizeEntryMatrix(bMatrix);
}

double DiisProcedure::bMatrixElement(int i, int j)
{
	int size = errorMatrix[0][0].size();
	double auxSum = 0.0e0;
	for (int ii = 0; ii < size; ii++)
	{
		for (int k = 0; k < size; k++)
		{
			auxSum += errorMatrix[i][ii][k] * errorMatrix[j][ii][k];
		}
	}
	return auxSum;
}


vector< vector<double> > DiisProcedure::matrix_product(const vector< vector<double> > &matrix_A, const vector< vector<double> > &matrix_B)
{
	vector< vector<double> > vec_out;
	int size = matrix_A.size();
	vec_out.resize(size);
	for (int is = 0; is<size; is++)
	{
		vec_out[is].resize(size);
	}

	for (int i_A = 0; i_A < size; i_A++)
	{
		for (int j_B = 0; j_B<size; j_B++)
		{
			vec_out[i_A][j_B] = 0.0e0;
			for (int k = 0; k < size; k++)
			{
				vec_out[i_A][j_B] += matrix_A[i_A][k] * matrix_B[k][j_B];
			}

		}
	}

	return vec_out;
}
