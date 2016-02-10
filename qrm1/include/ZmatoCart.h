/*
Artigo de referência:

Practical conversion from Torsion
Space to Cartesian Space for In
Silico Protein Synthesis
Jerod Parsons et.al.
DOI 10.1002/jcc.20237
*/

#ifndef ZMATOCART_H
#define ZMATOCART_H

#include <string>
#include <vector>

#include "Matrixop.h"
#include "Coordstructs.h"

class ZmatoCart: public Matrixop
{
	void ztocart1atomo(std::vector<CoordZMAT>);
	void ztocart2atomos(std::vector<CoordZMAT>);
	void ztocart3atomos(std::vector<CoordZMAT>);
	void ztocartiatomo(std::vector<CoordZMAT>,int);
public:
	std::vector<CoordXYZ> MolXYZ;
	std::vector<CoordZMAT> MolZMAT;
	void ztocart(std::vector<CoordZMAT>);
//	void ZmatoCart::operator=(std::vector<CoordZMAT>);
//	void ZmatoCart::operator=(std::vector<CoordXYZ>);
};

#endif
