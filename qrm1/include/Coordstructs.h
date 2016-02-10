#ifndef COORDSTRUCTS_H
#define COORDSTRUCTS_H

#include <string>

struct CoordXYZ
{
	std::string atomlabel;
	double x, y, z;
};

struct CoordZMAT
{
	std::string atomlabel;
	double dis;
	double ang;
	double die;
	int conect[3];
	int flagconec[3];
};

struct diatomicFourCenter
{
	std::vector< std::vector< std::vector< std::vector<double> > > > f_center;
};


#endif