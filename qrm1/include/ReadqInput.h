#ifndef READQINPUT_H
#define READQINPUT_H

#include <string>
#include <vector>

class ReadqInput
{
public:
	ReadqInput(std::string inputName);
	~ReadqInput();

	std::string getScfType();
	inline int getDebugLevel(){ return debugLevel; }
	std::string printLogName();
	inline int getCharge(){ return charge; }
	inline std::vector<double> getCoordinates() { return coordinates; }
	inline std::vector<std::string> getLabels(){ return labels; }
	inline int getParametricType(){ return parametricType; }
	inline bool getPrintExcell(){ return printExcell; }

private:
	int parametricType;
	int nAtoms;
	int charge;
	int scfType;
	int debugLevel;
	bool printExcell;
	std::vector<std::string> labels;
	std::vector<double> coordinates;
};

#endif
