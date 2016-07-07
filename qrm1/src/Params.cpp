#include "Params.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>

using namespace std;

void Params::set_semiempirical_parameters(string method_in)
{
	int number_of_registered_atoms = 3;
	int number_of_registered_double_parameters = 36;
	int number_of_registered_int_parameters = 4;

	set_parameter_vector_sizes(number_of_registered_atoms, number_of_registered_double_parameters, number_of_registered_int_parameters);
	hartree_ev = 27.2113961e0;
	ev_hartree = 1 / hartree_ev;
	bohr_angs = 0.529177249e0;
	angs_bohr = 1 / bohr_angs;
	pi = 3.14159265358979;

	if (method_in == "RM1")
	{
		method = "RM1";
		set_rm1_H_parameters();
		set_rm1_F_parameters();
		set_rm1_O_parameters();
	}
	else
	{
		cout << "METODO NAO ENCONTRADO" << endl
			<< "Nao foi possivel iniciar os parametros em Semiempirical_Parameters"
			<< endl;
		exit(3);
	}
}

double Params::set_d1_parameter(int principal_quantum_number, double expoents, double expoentp)
{
	double n_quantum_number = (double)principal_quantum_number;
	return ((1.0e0 + 2.0e0*n_quantum_number) / sqrt(3.0e0))*
		(
		pow(4.0e0*expoents*expoentp,
		n_quantum_number + 0.5e0) /
		pow(expoents + expoentp,
		2.0e0*n_quantum_number + 2.0e0)
		);
}

double Params::set_d2_parameter(int principal_quantum_number, double expoentp)
{
	double n_quantum_number = (double)principal_quantum_number;
	return sqrt(((1.0e0 + 2.0e0*n_quantum_number)*(2.0e0 + 2.0e0*n_quantum_number)) / 20.0e0) / expoentp;
}

double Params::set_additive_term_ro_0(double gss)
{
	return (1.0e0 / (2.0e0 * gss));
}


double Params::set_additive_term_ro_1(double hsp, double D1)
{
	/*
	Comes from an iterative method. This version was taken at mopac manual:http://openmopac.net/manual/
	Variable names are in agreement with the site. AD and a.
	*/

	double iter_ad_n2 = pow((D1*D1) / hsp, (-1.0e0 / 3.0e0));
	double iter_a_n2 = 0.5e0 * iter_ad_n2 - 0.5e0 *(1.0e0 /
		sqrt(4.0e0 * D1*D1 + 1.0e0 / (iter_ad_n2*iter_ad_n2)));
	double iter_ad_n1 = iter_ad_n2*(hsp / iter_a_n2);

	double iter_a_n1 = 0.0e0;
	double iter_ad_n = 0.0e0;
	double precisao_ad = 1.0e0;
	double precisao_a = 1.0e0;
	int i = 0;

	while (abs(precisao_ad) > 1.0e-8)
	{
		if (i != 0)
		{
			iter_a_n2 = iter_a_n1;
			iter_ad_n2 = iter_ad_n1;
			iter_ad_n1 = iter_ad_n;
		}

		iter_a_n1 = 0.5e0 * iter_ad_n1 - 0.5e0 *(1.0e0 /
			sqrt(4.0e0 * D1*D1 + 1.0e0 / (iter_ad_n1*iter_ad_n1)));

		precisao_ad = iter_ad_n1 - iter_ad_n2;
		precisao_a = iter_a_n1 - iter_a_n2;
		if (abs(precisao_a) < 1.0e-8)
		{
			break;
		}

		iter_ad_n = iter_ad_n2 + precisao_ad*
			((hsp - iter_a_n2) / (precisao_a));

		i++;

		if (i > 30)
		{
			cout << "ERRO NO ITERADOR DA INTEGRAL DE 4 CENTROS" << endl;
			exit(6);
		}
	}


	return 1.0e0 / (iter_ad_n*2.0e0);

	/*
	TEST TO VERIFY IF THIS CALCULATED RO IS CORRECT:
	float ev_ua = (float) 0.0367493088244753e0;
	float testa_ro_1_zz = 0.25e0*(1.0e0 / ro_1_B) - 0.5e0*(1.0e0 / sqrt(4.0e0*D1_B*D1_B + 4.0e0*ro_1_B*ro_1_B));
	cout << "testa ro1:  " << testa_ro_1_zz / ev_ua << endl;
	cout << "hsp:  " << mol.atomos_mol[B].bases_do_atomo[lambda].hsp / ev_ua << endl << endl;
	*/
}

double Params::set_additive_term_ro_2(double hpp, double D2)
{
	/*
	Comes from an iterative method. This version was taken at mopac manual:http://openmopac.net/manual/
	Variable names are in agreement with the site. AD and a.
	*/

	double iter_aq_n2 = pow(hpp / (3.0e0*D2*D2*D2*D2), (1.0e0 / 5.0e0));
	double iter_a_n2 = 0.25e0 * iter_aq_n2 - 0.5e0 *(1.0e0 /
		sqrt(4.0e0 * D2*D2 + 1.0e0 / (iter_aq_n2*iter_aq_n2))) +
		0.25e0*(1.0e0 / sqrt(8.0e0*D2*D2 + 1.0e0 / (iter_aq_n2*iter_aq_n2)));
	double iter_aq_n1 = iter_aq_n2*(hpp / iter_a_n2);

	double iter_a_n1 = 0.0e0;
	double iter_aq_n = 0.0e0;
	double precisao_aq = 1.0e0;
	double precisao_a = 1.0e0;
	int i = 0;

	while (abs(precisao_aq) > 1.0e-8)
	{
		if (i != 0)
		{
			iter_a_n2 = iter_a_n1;
			iter_aq_n2 = iter_aq_n1;
			iter_aq_n1 = iter_aq_n;
		}

		iter_a_n1 = 0.25e0 * iter_aq_n1 - 0.5e0 *(1.0e0 /
			sqrt(4.0e0 * D2*D2 + 1.0e0 / (iter_aq_n1*iter_aq_n1))) +
			0.25e0*(1.0e0 / sqrt(8.0e0*D2*D2 + 1.0e0 / (iter_aq_n1*iter_aq_n1)));

		precisao_aq = iter_aq_n1 - iter_aq_n2;
		precisao_a = iter_a_n1 - iter_a_n2;
		if (abs(precisao_a) < 1.0e-8)
		{
			break;
		}

		iter_aq_n = iter_aq_n2 + precisao_aq*
			((hpp - iter_a_n2) / (precisao_a));

		i++;

		if (i > 200)
		{
			cout << "ERRO NO ITERADOR DA INTEGRAL DE 4 CENTROS" << endl;
			exit(6);
		}
	}

	// test to verify if it is convergint to the right point
	/*
	double ro_2 = 1.0e0 / (iter_aq_n*2.0e0);
	double test = 0.25e0*pow(4.0e0*ro_2*ro_2, -0.5e0)
	+ 0.25e0*pow(8.0e0*D2*D2 + 4.0e0*ro_2*ro_2, -0.5e0)
	- 0.5e0*pow(4.0e0*D2*D2 + 4.0e0*ro_2*ro_2, -0.5e0);
	cout << "valor correto:  " << hpp << "  valor do teste:  " << test << endl;
	*/

	return 1.0e0 / (iter_aq_n*2.0e0);

}

/////////////////////////
// MANAGING PARAMETERS //
/////////////////////////
double Params::get_double(string atom_type, string parameter)
{

	if (atom_type == "H")
	{
		return double_parameters[0][get_number_parameter_double(parameter)];
	}
	else if (atom_type == "F")
	{
		return double_parameters[1][get_number_parameter_double(parameter)];
	}
	else if(atom_type == "O")
	{
		return double_parameters[2][get_number_parameter_double(parameter)];
	}
	else
	{
		cout << "atomo ainda nao cadastrado - erro em Semiempirical_Parameters" << endl;
		exit(3);
	}

}

int Params::get_number_parameter_double(string parameter)
{
	if (parameter == "uss")
	{
		return 0;
	}
	else if (parameter == "upp")
	{
		return 1;
	}
	else if (parameter == "betas")
	{
		return 2;
	}
	else if (parameter == "betap")
	{
		return 3;
	}
	else if (parameter == "alfacore")
	{
		return 4;
	}
	else if (parameter == "gss")
	{
		return 5;
	}
	else if (parameter == "gsp")
	{
		return 6;
	}
	else if (parameter == "gpp")
	{
		return 7;
	}
	else if (parameter == "gp2")
	{
		return 8;
	}
	else if (parameter == "hsp")
	{
		return 9;
	}
	else if (parameter == "a1core")
	{
		return 10;
	}
	else if (parameter == "b1core")
	{
		return 11;
	}
	else if (parameter == "c1core")
	{
		return 12;
	}
	else if (parameter == "a2core")
	{
		return 13;
	}
	else if (parameter == "b2core")
	{
		return 14;
	}
	else if (parameter == "c2core")
	{
		return 15;
	}
	else if (parameter == "a3core")
	{
		return 16;
	}
	else if (parameter == "b3core")
	{
		return 17;
	}
	else if (parameter == "c3core")
	{
		return 18;
	}
	else if (parameter == "a4core")
	{
		return 19;
	}
	else if (parameter == "b4core")
	{
		return 20;
	}
	else if (parameter == "c4core")
	{
		return 21;
	}
	else if (parameter == "expoents")
	{
		return 22;
	}
	else if (parameter == "expoentp")
	{
		return 23;
	}
	else if (parameter == "D1")
	{
		return 24;
	}
	else if (parameter == "D2")
	{
		return 25;
	}
	else if (parameter == "ro_0")
	{
		return 26;
	}
	else if (parameter == "ro_1")
	{
		return 27;
	}
	else if (parameter == "ro_2")
	{
		return 28;
	}
	else if (parameter == "qover")
	{
		return 29;
	}
	else if (parameter == "qmono")
	{
		return 30;
	}
	else if (parameter == "qalfa")
	{
		return 31;
	}
	else if (parameter == "qgauss")
	{
		return 32;
	}
	else if (parameter == "qr0")
	{
		return 33;
	}
	else if (parameter == "qrinf")
	{
		return 34;
	}
	else if (parameter == "skew")
	{
		return 35;
	}
	else
	{
		cout << "parametro nao cadastrado" << endl
			<< " verificar a rotina Semiempirical_Parameters" << endl;
		exit(3);
	}
}

int Params::get_int(string atom_type, string parameter)
{

	if (atom_type == "H")
	{
		return int_parameters[0][get_number_parameter_int(parameter)];
	}
	else if (atom_type == "F")
	{
		return int_parameters[1][get_number_parameter_int(parameter)];
	}
	else if (atom_type == "O")
	{
		return int_parameters[2][get_number_parameter_int(parameter)];
	}
	else
	{
		cout << "atomo ainda nao cadastrado - erro em Semiempirical_Parameters" << endl;
		exit(3);
	}

}

int Params::get_number_parameter_int(string parameter)
{
	if (parameter == "q_number")
	{
		return 0;
	}
	else if (parameter == "base_number")
	{
		return 1;
	}
	else if (parameter == "number_of_electrons")
	{
		return 2;
	}
	else if (parameter == "atomic_charge")
	{
		return 3;
	}
	else
	{
		cout << "parametro nao cadastrado - consultar o objeto Semiempirical_Parameters"
			<< endl;
		exit(3);
	}
}


void Params::set_parameter_vector_sizes(int number_of_registered_atoms, int number_of_registered_double_parameters, int number_of_registered_int_parameters)
{
	double_parameters.resize(number_of_registered_atoms);
	int_parameters.resize(number_of_registered_atoms);
	for (int i = 0; i < number_of_registered_atoms; i++)
	{
		double_parameters[i].resize(number_of_registered_double_parameters);
		int_parameters[i].resize(number_of_registered_int_parameters);
	}
}

void Params::set_rm1_H_parameters()
{
	// parametros eletronicos em ua
	// parametros nucleares em eV
	int_parameters[0][0] = 1; //principal_quantum_number
	int_parameters[0][1] = 1; //base_number
	int_parameters[0][2] = 1; //number_of_electrons
	int_parameters[0][3] = 1; // atomic_number

	double_parameters[0][0] = -11.96067700e0*ev_hartree; //uss
	double_parameters[0][1] = 0.0e0; //upp
	double_parameters[0][2] = -5.76544470e0*ev_hartree; //betas 
	double_parameters[0][3] = 0.0e0; //betap
	double_parameters[0][4] = 3.06835950e0; //alfacore
	double_parameters[0][5] = 13.98321300e0*ev_hartree; //gss
	double_parameters[0][6] = 0.0e0; //gsp
	double_parameters[0][7] = 0.0e0; //gpp
	double_parameters[0][8] = 0.0e0; //gp2
	double_parameters[0][9] = 0.0e0; //hsp
	double_parameters[0][10] = 0.10288880e0; //a1core
	double_parameters[0][11] = 5.90172270e0; //b1core
	double_parameters[0][12] = 1.17501180e0; //c1core
	double_parameters[0][13] = 0.06457450e0; //a2core
	double_parameters[0][14] = 6.41785670e0; //b2core
	double_parameters[0][15] = 1.93844480e0; //c2core
	double_parameters[0][16] = -0.03567390e0; //a3core
	double_parameters[0][17] = 2.80473130e0; //b3core
	double_parameters[0][18] = 1.63655240e0; //c3core
	double_parameters[0][19] = 0.0e0; //a4core
	double_parameters[0][20] = 0.0e0; //b4core
	double_parameters[0][21] = 0.0e0; //c4core
	double_parameters[0][22] = 1.08267370e0; //expoents
	double_parameters[0][23] = 0.0e0; //expoentp
	double_parameters[0][24] = 0.0e0; //D1
	double_parameters[0][25] = 0.0e0; //D2
	double_parameters[0][26] = set_additive_term_ro_0(double_parameters[0][5]); //ro_0
	double_parameters[0][27] = 0.0e0; //ro_1
	double_parameters[0][28] = 0.0e0; //ro_2
	double_parameters[0][29] = 1.0e0; // q - overlap
	double_parameters[0][30] = 1.0e0; // q - mono
	double_parameters[0][31] = 1.3e0; // q - alfa
	double_parameters[0][32] = 1.0e0; // q - gauss
	double_parameters[0][33] = 0.0e0; // integral - q - R==0
	double_parameters[0][34] = 0.0e0; // integral - q - R==inf
	double_parameters[0][35] = 0.0e0; // gaussian skew

}

void Params::set_qrm1_H_parameters()
{
	// parametros eletronicos em ua
	// parametros nucleares em eV
	int_parameters[0][0] = 1; //principal_quantum_number
	int_parameters[0][1] = 1; //base_number
	int_parameters[0][2] = 1; //number_of_electrons
	int_parameters[0][3] = 1; // atomic_number

	double_parameters[0][0] = -11.96067700e0*ev_hartree; //uss
	double_parameters[0][1] = 0.0e0; //upp
	double_parameters[0][2] = -5.76544470e0*ev_hartree; //betas 
	double_parameters[0][3] = 0.0e0; //betap
	double_parameters[0][4] = 3.06835950e0; //alfacore
	double_parameters[0][5] = 13.98321300e0*ev_hartree; //gss
	double_parameters[0][6] = 0.0e0; //gsp
	double_parameters[0][7] = 0.0e0; //gpp
	double_parameters[0][8] = 0.0e0; //gp2
	double_parameters[0][9] = 0.0e0; //hsp
	double_parameters[0][10] = 0.10288880e0; //a1core
	double_parameters[0][11] = 5.90172270e0; //b1core
	double_parameters[0][12] = 1.17501180e0; //c1core
	double_parameters[0][13] = 0.06457450e0; //a2core
	double_parameters[0][14] = 6.41785670e0; //b2core
	double_parameters[0][15] = 1.93844480e0; //c2core
	double_parameters[0][16] = -0.03567390e0; //a3core
	double_parameters[0][17] = 2.80473130e0; //b3core
	double_parameters[0][18] = 1.63655240e0; //c3core
	double_parameters[0][19] = 0.0e0; //a4core
	double_parameters[0][20] = 0.0e0; //b4core
	double_parameters[0][21] = 0.0e0; //c4core
	double_parameters[0][22] = 1.08267370e0; //expoents
	double_parameters[0][23] = 0.0e0; //expoentp
	double_parameters[0][24] = 0.0e0; //D1
	double_parameters[0][25] = 0.0e0; //D2
	double_parameters[0][26] = set_additive_term_ro_0(double_parameters[0][5]); //ro_0
	double_parameters[0][27] = 0.0e0; //ro_1
	double_parameters[0][28] = 0.0e0; //ro_2
	double_parameters[0][29] = 1.001e0; // q - overlap
	double_parameters[0][30] = 1.001e0; // q - monopole
}


void Params::set_rm1_F_parameters()
{
	// parametros eletronicos em ua
	// parametros nucleares em eV
	int_parameters[1][0] = 2; //principal_quantum_number
	int_parameters[1][1] = 4; //base_number
	int_parameters[1][2] = 7; //number_of_electrons
	int_parameters[1][3] = 7; // atomic_charge

	double_parameters[1][0] = -134.18369591e0*ev_hartree; //uss
	double_parameters[1][1] = -107.84660920e0*ev_hartree; //upp
	double_parameters[1][2] = -70.00000512e0*ev_hartree; //betas 
	double_parameters[1][3] = -32.67982711e0*ev_hartree; //betap
	double_parameters[1][4] = 6.00000062e0; //alfacore
	double_parameters[1][5] = 16.72091319e0*ev_hartree; //gss
	double_parameters[1][6] = 16.76142629e0*ev_hartree; //gsp
	double_parameters[1][7] = 15.22581028e0*ev_hartree; //gpp
	double_parameters[1][8] = 14.86578679e0*ev_hartree; //gp2
	double_parameters[1][9] = 1.99766171e0*ev_hartree; //hsp
	double_parameters[1][10] = 0.40302025e0; //a1core
	double_parameters[1][11] = 7.20441959e0; //b1core
	double_parameters[1][12] = 0.81653013e0; //c1core
	double_parameters[1][13] = 0.07085831e0; //a2core
	double_parameters[1][14] = 9.00001562e0; //b2core
	double_parameters[1][15] = 1.43802381e0; //c2core
	double_parameters[1][16] = 0.0e0; //a3core
	double_parameters[1][17] = 0.0e0; //b3core
	double_parameters[1][18] = 0.0e0; //c3core
	double_parameters[1][19] = 0.0e0; //a4core
	double_parameters[1][20] = 0.0e0; //b4core
	double_parameters[1][21] = 0.0e0; //c4core
	double_parameters[1][22] = 4.40337913e0; //expoents
	double_parameters[1][23] = 2.64841556e0; //expoentp
	double_parameters[1][24] = set_d1_parameter(int_parameters[1][0], double_parameters[1][22], double_parameters[1][23]); //D1
	double_parameters[1][25] = set_d2_parameter(int_parameters[1][0], double_parameters[1][23]); //D2
	double_parameters[1][26] = set_additive_term_ro_0(double_parameters[1][5]); //ro_0
	double_parameters[1][27] = set_additive_term_ro_1(double_parameters[1][9], double_parameters[1][24]); //ro_1
	double_parameters[1][28] = set_additive_term_ro_2(0.5e0*(double_parameters[1][7] - double_parameters[1][8]), double_parameters[1][25]); //ro_2
}


void Params::set_rm1_O_parameters()
{
	// parametros eletronicos em ua
	// parametros nucleares em eV
	int_parameters[2][0] = 2; //principal_quantum_number
	int_parameters[2][1] = 4; //base_number
	int_parameters[2][2] = 6; //number_of_electrons
	int_parameters[2][3] = 6; // atomic_charge

	double_parameters[2][0] = -96.94948070e0*ev_hartree; //uss
	double_parameters[2][1] = -77.89092980e0*ev_hartree; //upp
	double_parameters[2][2] = -29.85101210e0*ev_hartree; //betas 
	double_parameters[2][3] = -29.15101310e0*ev_hartree; //betap
	double_parameters[2][4] = 4.17196720e0; //alfacore
	double_parameters[2][5] = 14.00242790e0*ev_hartree; //gss
	double_parameters[2][6] = 14.95625040e0*ev_hartree; //gsp
	double_parameters[2][7] = 14.14515140e0*ev_hartree; //gpp
	double_parameters[2][8] = 12.70325500e0*ev_hartree; //gp2
	double_parameters[2][9] = 3.93217160e0*ev_hartree; //hsp
	double_parameters[2][10] = 0.23093550e0; //a1core
	double_parameters[2][11] = 5.21828740e0; //b1core
	double_parameters[2][12] = 0.90363560e0; //c1core
	double_parameters[2][13] = 0.05859870e0; //a2core
	double_parameters[2][14] = 7.42932930e0; //b2core
	double_parameters[2][15] = 1.51754610e0; //c2core
	double_parameters[2][16] = 0.0e0; //a3core
	double_parameters[2][17] = 0.0e0; //b3core
	double_parameters[2][18] = 0.0e0; //c3core
	double_parameters[2][19] = 0.0e0; //a4core
	double_parameters[2][20] = 0.0e0; //b4core
	double_parameters[2][21] = 0.0e0; //c4core
	double_parameters[2][22] = 3.17936910e0; //expoents
	double_parameters[2][23] = 2.55361910e0; //expoentp
	double_parameters[2][24] = set_d1_parameter(int_parameters[2][0], double_parameters[2][22], double_parameters[2][23]); //D1
	double_parameters[2][25] = set_d2_parameter(int_parameters[2][0], double_parameters[2][23]); //D2
	double_parameters[2][26] = set_additive_term_ro_0(double_parameters[2][5]); //ro_0
	double_parameters[2][27] = set_additive_term_ro_1(double_parameters[2][9], double_parameters[2][24]); //ro_1
	double_parameters[2][28] = set_additive_term_ro_2(0.5e0*(double_parameters[2][7] - double_parameters[2][8]), double_parameters[2][25]); //ro_2


//	cout.precision(8);
//	cout << " d1:  " << double_parameters[2][24] << endl;
//	cout << " d2:  " << double_parameters[2][25] << endl;
//	cout << " ro0:  " << double_parameters[2][26] << endl;
//	cout << " ro1:  " << double_parameters[2][27] << endl;
//	cout << " ro2:  " << double_parameters[2][28] << endl;
//	cin.get();

}



// AAAAMMMMM11111
void Params::set_am1_O_parameters()
{
	// parametros eletronicos em ua
	// parametros nucleares em eV
	int_parameters[2][0] = 2; //principal_quantum_number
	int_parameters[2][1] = 4; //base_number
	int_parameters[2][2] = 6; //number_of_electrons
	int_parameters[2][3] = 6; // atomic_charge

	double_parameters[2][0] = -97.8300000e0*ev_hartree; //uss
	double_parameters[2][1] = -78.26238000e0*ev_hartree; //upp
	double_parameters[2][2] = -29.27277300e0*ev_hartree; //betas 
	double_parameters[2][3] = -29.27277300e0*ev_hartree; //betap
	double_parameters[2][4] = 4.17196717e0; //alfacore
	double_parameters[2][5] = 15.42000000e0*ev_hartree; //gss
	double_parameters[2][6] = 14.48000000e0*ev_hartree; //gsp
	double_parameters[2][7] = 14.52000000e0*ev_hartree; //gpp
	double_parameters[2][8] = 12.98000000e0*ev_hartree; //gp2
	double_parameters[2][9] = 3.94000000e0*ev_hartree; //hsp
	double_parameters[2][10] = 0.23093552e0; //a1core
	double_parameters[2][11] = 5.21828736e0; //b1core
	double_parameters[2][12] = 0.90363555e0; //c1core
	double_parameters[2][13] = 0.05859873e0; //a2core
	double_parameters[2][14] = 7.42932932e0; //b2core
	double_parameters[2][15] = 1.51754610e0; //c2core
	double_parameters[2][16] = 0.0e0; //a3core
	double_parameters[2][17] = 0.0e0; //b3core
	double_parameters[2][18] = 0.0e0; //c3core
	double_parameters[2][19] = 0.0e0; //a4core
	double_parameters[2][20] = 0.0e0; //b4core
	double_parameters[2][21] = 0.0e0; //c4core
	double_parameters[2][22] = 3.10803200e0; //expoents
	double_parameters[2][23] = 2.52403900e0; //expoentp
	double_parameters[2][24] = set_d1_parameter(int_parameters[2][0], double_parameters[2][22], double_parameters[2][23]); //D1
	double_parameters[2][25] = set_d2_parameter(int_parameters[2][0], double_parameters[2][23]); //D2
	double_parameters[2][26] = set_additive_term_ro_0(double_parameters[2][5]); //ro_0
	double_parameters[2][27] = set_additive_term_ro_1(double_parameters[2][9], double_parameters[2][24]); //ro_1
	double_parameters[2][28] = set_additive_term_ro_2(0.5e0*(double_parameters[2][7] - double_parameters[2][8]), double_parameters[2][25]); //ro_2
}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////// PARAMETRIZACAO  ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


void Params::readFromFile(string fileName)
{
	cout << "READ FROM FILE PRECISA DE REFORMAS" << endl;
	exit(3);
}

void Params::readFromFileParametrization(string fileName)
{
	ifstream paramFile_;
	paramFile_.open(fileName.c_str());

	int methodNumber;
	paramFile_ >> methodNumber;
	method = getMethod(methodNumber);

	readFileMethod(paramFile_);
	paramFile_.close();
}

string Params::getMethod(int methodNumber)
{
	switch (methodNumber)
	{
	case 0:
		return "RM1-3g";
		break;

	case 1:
		return "RM1-2g";
		break;

	case 2:
		return "RM1-1g";
		break;

	case 3:
		return "qRM1-3g";
		break;

	case 4:
		return "qRM1-2g";
		break;

	case 5:
		return "qRM1-1g";
		break;

        case 6:
		return "skewRM1";
		break;

	default:
		cout << "Erro em getMethod" << endl;
		exit(3);
		break;
	}
}


void Params::readFileMethod(ifstream &paramFile_)
{
	double readParam;
	paramFile_ >> readParam;
//	double_parameters[0][0] = readParam*ev_hartree;
	double_parameters[0][0] *= readParam;
	paramFile_ >> readParam;
//	double_parameters[0][2] = readParam*ev_hartree;
	double_parameters[0][2] *= readParam;
	paramFile_ >> readParam;
//	double_parameters[0][4] = readParam;
	double_parameters[0][4] *= readParam;
	paramFile_ >> readParam;
//	double_parameters[0][5] = readParam*ev_hartree;
	double_parameters[0][5] *= readParam;

	// na nova versao zeta e aqui!!!
	//zeta
	paramFile_ >> readParam;
	//	double_parameters[0][22] = readParam;
	double_parameters[0][22] *= readParam;
	double_parameters[0][26] = set_additive_term_ro_0(double_parameters[0][5]); //ro_0



	paramFile_ >> readParam;
//	double_parameters[0][10] = readParam;
	double_parameters[0][10] *= readParam;
	paramFile_ >> readParam;
//	double_parameters[0][11] = readParam;
	double_parameters[0][11] *= readParam;
	paramFile_ >> readParam;
//	double_parameters[0][12] = readParam;
	double_parameters[0][12] *= readParam;

	if (method == "skewRM1")
	{
		paramFile_ >> readParam;
		double_parameters[0][13] = readParam;
		paramFile_ >> readParam;
		double_parameters[0][14] = readParam;
		paramFile_ >> readParam;
		double_parameters[0][15] = readParam;
		paramFile_ >> readParam;
		double_parameters[0][35] = readParam; //skew
	}
	

	// parando aqui e o RM1-1g
	if (
		(method == "RM1-2g") || (method == "RM1-3g")||
		(method == "qRM1-2g") || (method == "qRM1-3g")
		)
	{
		paramFile_ >> readParam;
//		double_parameters[0][13] = readParam;
		double_parameters[0][13] *= readParam;
		paramFile_ >> readParam;
//		double_parameters[0][14] = readParam;
		double_parameters[0][14] *= readParam;
		paramFile_ >> readParam;
//		double_parameters[0][15] = readParam;
		double_parameters[0][15] *= readParam;
	}
	if (
		(method == "RM1-3g") ||
		(method == "qRM1-3g")
		)
	{
		paramFile_ >> readParam;
//		double_parameters[0][16] = readParam;
		double_parameters[0][16] *= readParam;
		paramFile_ >> readParam;
//		double_parameters[0][17] = readParam;
		double_parameters[0][17] *= readParam;
		paramFile_ >> readParam;
//		double_parameters[0][18] = readParam;
		double_parameters[0][18] *= readParam;
	}
	if (
		(method == "qRM1-1g") ||
		(method == "qRM1-2g") ||
		(method == "qRM1-3g")
		)
	{
		paramFile_ >> readParam;
		double_parameters[0][29] = readParam;//qover
		paramFile_ >> readParam;
		double_parameters[0][30] = readParam;//qmono
		paramFile_ >> readParam;
		double_parameters[0][31] = readParam;//qalfa
		paramFile_ >> readParam;
		double_parameters[0][32] = readParam;//qgauss
	}

	// SETTING qr0 and qinf
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double ro_0 = double_parameters[0][26];

	double q = double_parameters[0][30];
	if (
		(method[0] == 'q') &&
		(abs(q - 1.0e0) > 0.000001)
		)
	{
		double x = 4.0e0*ro_0*ro_0;
		double_parameters[0][34] = pow(1.5e0 - 0.5e0*pow((x+1.0e12), 1.0e0 - q), (1.0e0 / (1.0e0 - q)));
		double_parameters[0][33] = pow(1.5e0 - 0.5e0*pow(x, 1.0e0 - q), (1.0e0 / (1.0e0 - q)))
			- double_parameters[0][34];
	}
	else
	{
		double_parameters[0][33] = double_parameters[0][5];
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// DONE


}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////// FIM DA PARAMETRIZACAO  ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////





// VARIAVEIS STATICAS PRECISAM DE DEFINICAO EM TODOS OS ARQUIVOS
double Params::ev_hartree;
double Params::hartree_ev;
double Params::bohr_angs;
double Params::angs_bohr;
double Params::pi;
string Params::method;

vector< vector<int> > Params::int_parameters;
vector< vector<double> > Params::double_parameters;
