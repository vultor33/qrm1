/*
void Ciclo_SCF::adiciona_shift_aleatorio()
{
	float Energia_eletronica = 0.0e0;
	float X = 2.0e0;
	//szabo pag 150
	for (int mi = 0; mi < fock_matrix_size; mi++)
	{
		for (int ni = 0; ni < fock_matrix_size; ni++)
		{
			float r1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / X));
			fock_matrix[mi][ni] += r1;
		}
	}

}

void Ciclo_SCF::adiciona_shift_b(float b)
{
	for (int mi = 0; mi < fock_matrix_size; mi++)
	{
		for (int ni = 0; ni < fock_matrix_size; ni++)
		{
			fock_matrix[mi][ni] += b*density_matrix[mi][ni];
		}
	}
}





// Violates rotational invariance
// Its integral 22 of thiel article. We resolve this with: 0,5*((yy|yy)-(yy|zz))
double Four_Center_Integrals::Qyz_Qyz_bracket(const Molecule &mol, int A, int B)
{
	double RAB = mol.distancia_i_j[A][B];

	double D2_A = PARAM::get_double(mol.atom_name[A], "D2");
	double ro_2_A = PARAM::get_double(mol.atom_name[A], "ro_2");
	double D2_B = PARAM::get_double(mol.atom_name[B], "D2");
	double ro_2_B = PARAM::get_double(mol.atom_name[B], "ro_2");

	double a_2_2_quad = (ro_2_A + ro_2_B)*(ro_2_A + ro_2_B);

	double colchete =  0.25e0*(1.0e0 / sqrt(RAB*RAB + 2.0e0*(D2_A - D2_B)*(D2_A - D2_B) + a_2_2_quad))
		+ 0.25e0*(1.0e0 / sqrt(RAB*RAB + 2.0e0*(D2_A + D2_B)*(D2_A + D2_B) + a_2_2_quad))
		- 0.5e0*(1.0e0 / sqrt(RAB*RAB + 2.0e0*D2_A*D2_A + 2.0e0*D2_B*D2_B + a_2_2_quad))
		;
	return colchete;
}




07/05/2015
void ScfCycle::build_fock_matrix()
{
int i_fock = -1;
int j_fock = -1;
int auxni;
for (int A = 0; A < mol.number_of_atoms; A++)
{
int n_bases_A = Params::get_int(mol.atom_name[A], "base_number");
for (int mi = 0; mi < n_bases_A; mi++)
{
i_fock++;
j_fock = i_fock;
for (int B = A; B < mol.number_of_atoms; B++)
{
int n_bases_B = Params::get_int(mol.atom_name[B], "base_number");
if (B == A)
{
auxni = mi;
}
else
{
auxni = 0;
}
for (int ni = auxni; ni < n_bases_B; ni++)
{
if ((A == B) && (mi == ni))
{
//densidade mi mi
fock_matrix[i_fock][j_fock] += couloumb_exchange_mi_mi(A, mi);
}
else if ((A == B))
{
fock_matrix[i_fock][j_fock] += couloumb_exchange_mi_ni(A, mi, ni);
}
else
{
fock_matrix[i_fock][j_fock] += couloumb_exchange_mi_lambda(A, B, mi, ni);
}
j_fock++;
}
}
}
}
}

double ScfCycle::couloumb_exchange_mi_lambda(int A, int B, int mi, int lambda)
{
// ni e o lambda nesse caso

double auxsoma = 0.0e0;
int i_densidade;
int j_densidade;

for (int ni = 0; ni < Params::get_int(mol.atom_name[A],"base_number"); ni++)
{
for (int sigma = 0; sigma < Params::get_int(mol.atom_name[B], "base_number"); sigma++)
{
i_densidade = mol.i_fock_base_line[A][ni];
j_densidade = mol.i_fock_base_line[B][sigma];
auxsoma += density_matrix[i_densidade][j_densidade] *
four_center_.get_four_center(A, B, mi, ni, lambda, sigma);

cout << "densi  " << density_matrix[i_densidade][j_densidade] << endl;
}
}

auxsoma *= -0.5e0;
return auxsoma;
}


double ScfCycle::couloumb_exchange_mi_ni(int A, int mi, int ni)
{
double auxsoma = 0.0e0;

int i_densidade = mol.i_fock_base_line[A][mi];
int j_densidade = mol.i_fock_base_line[A][ni];

auxsoma += 0.5e0*density_matrix[i_densidade][j_densidade] *
(3.0e0 * four_center_.get_four_center(A, A, mi, ni, mi, ni) -
four_center_.get_four_center(A, A, mi, mi, ni, ni));

for (int B = 0; B < mol.number_of_atoms; B++)
{
if (B != A)
{
for (int lambda = 0; lambda < Params::get_int(mol.atom_name[B],"base_number"); lambda++)
{
for (int sigma = 0; sigma < Params::get_int(mol.atom_name[B], "base_number"); sigma++)
{
int i_densidade = mol.i_fock_base_line[B][lambda];
int j_densidade = mol.i_fock_base_line[B][sigma];

auxsoma += density_matrix[i_densidade][j_densidade] *
four_center_.get_four_center(A, B, mi, ni, lambda, sigma);
}
}
}
}

return auxsoma;
}

double ScfCycle::couloumb_exchange_mi_mi(int A, int mi)
{
double auxsoma = 0;
int i_densidade;
int j_densidade;

for (int ni = 0; ni < Params::get_int(mol.atom_name[A],"base_number"); ni++)
{
auxsoma += density_matrix[mol.i_fock_base_line[A][ni]][mol.i_fock_base_line[A][ni]] *
(four_center_.get_four_center(A, A, mi, mi, ni, ni) -
0.5e0*four_center_.get_four_center(A, A, mi, ni, mi, ni));
}


for (int B = 0; B < mol.number_of_atoms; B++)
{
if (B != A)
{
for (int lambda = 0; lambda < Params::get_int(mol.atom_name[B],"base_number"); lambda++)
{
for (int sigma = 0; sigma < Params::get_int(mol.atom_name[B], "base_number"); sigma++)
{
i_densidade = mol.i_fock_base_line[B][lambda];
j_densidade = mol.i_fock_base_line[B][sigma];
auxsoma += density_matrix[i_densidade][j_densidade] *
four_center_.get_four_center(A, B, mi, mi, lambda, sigma);

cout << "lam:  " << lambda << "sigma " << sigma << "  int  " << four_center_.get_four_center(A, B, mi, mi, lambda, sigma)*Params::hartree_ev << endl;

				}
			}
		}
	}

	return auxsoma;
}




TESTE DA MONTAGEM DA MATRIZ DE FOCK DO MEU

density_matrix[0][0] = 0.959364;
density_matrix[0][1] = 0.684932;
density_matrix[0][2] = 0.154701;
density_matrix[0][3] = 0.582608;
density_matrix[0][4] = 0.362319;
density_matrix[0][5] = 0.185947;
density_matrix[1][1] = 0.918870;
density_matrix[1][2] = 0.486070;
density_matrix[1][3] = 0.303752;
density_matrix[1][4] = -0.297901;
density_matrix[1][5] = -0.327125;
density_matrix[2][2] = 1.401301;
density_matrix[2][3] = -0.670924;
density_matrix[2][4] = -0.159004;
density_matrix[2][5] = 0.321464;
density_matrix[3][3] = 0.807339;
density_matrix[3][4] = 0.194216;
density_matrix[3][5] = -0.208148;
density_matrix[4][4] = 0.926465;
density_matrix[4][5] = 0.843546;
density_matrix[5][5] = 0.986660;

int n_bases = fock_matrix_size;
for (int i = 0; i < n_bases; i++)
{
for (int j = 0; j < n_bases; j++)
{
if (j < i)
{
density_matrix[i][j] = density_matrix[j][i];
}
}
}



for (int i = 0; i < number_of_atoms; i++)
{
for (int j = 0; j < number_of_atoms; j++)
{
if (i != j)
{
distancia_i_j[i][j] = Params::angs_bohr*
sqrt((coord[i] - coord[j])*(coord[i] - coord[j])

+ (coord[i + number_of_atoms] - coord[j + number_of_atoms])*
(coord[i + number_of_atoms] - coord[j + number_of_atoms])

+ (coord[i + 2 * number_of_atoms] - coord[j + 2 * number_of_atoms])*
(coord[i + 2 * number_of_atoms] - coord[j + 2 * number_of_atoms]));
}
else
{
distancia_i_j[i][j] = 0.0e0;
}
}
}






vector<double> xyz(12);
xyz[0] = 0.0e0;
xyz[1] = 1.4827e0;
xyz[2] = 0.3795e0;
xyz[3] = -0.3763e0;
xyz[4] = 0.0e0;
xyz[5] = 0.0e0;
xyz[6] = 2.4538e0;
xyz[7] = -1.0337e0;
xyz[8] = 0.0e0;
xyz[9] = 0.0e0;
xyz[10] = 0.0e0;
xyz[11] = 1.0565e0;

*/



/*
ifstream paramFile_;
paramFile_.open(fileName.c_str());

string type;
paramFile_ >> type;
if (type == "H")
method = "RM1";
else if (type == "qH")
method = "qRM1";

for (int i = 0; i < 31; i++)
{
if ((method == "RM1") && (i == 29))
break;

paramFile_ >> double_parameters[0][i];
}

double_parameters[0][0] *= ev_hartree;
double_parameters[0][2] *= ev_hartree;
double_parameters[0][5] *= ev_hartree;
double_parameters[0][26] = set_additive_term_ro_0(double_parameters[0][5]);

paramFile_.close();





void save_calculated_integral(double integral_AB, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);

void FourCenterIntegrals::save_calculated_integral(double integral_AB, int atomo_A, int atomo_B, int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma)
{
vector<int> marcador(6);
marcador[0] = atomo_A;
marcador[1] = atomo_B;
marcador[2] = orbital_mi;
marcador[3] = orbital_ni;
marcador[4] = orbital_lambda;
marcador[5] = orbital_sigma;
flag_of_calculated_integrals.push_back(marcador);
calculated_integrals.push_back(integral_AB);
}










*/





