#ifndef DLIBLINEARSYSTEM_H
#define DLIBLINEARSYSTEM_H

#include <vector>

class DlibLinearSystem
{
public:
	DlibLinearSystem(){}
	~DlibLinearSystem();

//	void teste();
	bool calculateBCoefficients(std::vector< std::vector<double> > &bMatrix);
	std::vector<double> getCoefficientsOut();

private:
	std::vector<double> getbVector(int size);
	std::vector<double> coefficientsOut;

};

#endif