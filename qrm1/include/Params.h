#ifndef Params_H
#define Params_H

#include <string>
#include <vector>
#include <fstream>

class Params
{
private:
	static std::vector< std::vector<int> > int_parameters;
	static std::vector< std::vector<double> > double_parameters;
	static void set_parameter_vector_sizes(int number_of_registered_atoms, int number_of_registered_double_parameters, int number_of_registered_int_parameters);
	static int get_number_parameter_double(std::string parameter);
	static int get_number_parameter_int(std::string parameter);

	static double set_d1_parameter(int principal_quantum_number, double expoents, double expoentp);
	static double set_d2_parameter(int principal_quantum_number, double expoentp);
	static double set_additive_term_ro_0(double gss);
	static double set_additive_term_ro_1(double hsp, double D1);
	static double set_additive_term_ro_2(double hpp, double D2);

	static void set_rm1_H_parameters();
	static void set_rm1_F_parameters();
	static void set_rm1_O_parameters();
	static void set_qrm1_H_parameters();
	static void set_am1_O_parameters();

	static void readFileMethod(std::ifstream &paramFile_);
	static std::string getMethod(int methodNumber);

public:
	static double ev_hartree;
	static double hartree_ev;
	static double bohr_angs;
	static double angs_bohr;
	static double pi;
	static std::string method;

	static void set_semiempirical_parameters(std::string method_in);
	static void readFromFile(std::string fileName);
	static void readFromFileParametrization(std::string fileName);
	static double get_double(std::string atom_type, std::string parameter);
	static int get_int(std::string atom_type, std::string parameter);

};

#endif


