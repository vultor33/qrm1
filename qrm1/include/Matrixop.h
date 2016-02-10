#ifndef MatrixopH
#define MatrixopH

#include <vector>

#include "Coordstructs.h"

class Matrixop
{
public:
	double produtescalar(double, double, double, double, double, double);
	std::vector<double> produtvetor(double, double, double, double, double, double);
	std::vector<double> normalizar(double, double, double);
	std::vector<double> somavec(double, double, double, double, double, double);
	std::vector<double> vecAB(double, double, double, double, double, double);
	std::vector<double> matrizvec(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>);
	std::vector< std::vector<double> > matrix_matrix(std::vector< std::vector<double> > matA, std::vector< std::vector<double> > matB);// square matrix
};


#endif


//Rotina com varias operacoes entre vetores e matrizes uteis
