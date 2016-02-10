#ifndef CHECKTYPEOFINTEGRALTOBECALCULATED_H
#define CHECKTYPEOFINTEGRALTOBECALCULATED_H

class CheckTypeOfIntegralToBeCalculated
{
private:
	//one center conditions
	// so tem dois iguais (mas preciso colocar um exit para garantir)
	bool is_p_orbital(int orbital); //true if orb == p
	bool ss_ss_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool ss_pp_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);// spsp or pssp or spps or psps
	bool sp_sp_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);// spsp or pssp or spps or psps
	bool pp_pp_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool pp_plpl_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool ppl_ppl_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);

	// rest of two center conditions
	bool ss_spz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);//spx or pxs
	bool spz_ss_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);//spx or pxs
	bool ss_pzpz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool pzpz_ss_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool ss_ppippi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool ppippi_ss_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool sppi_sppi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);// (sppi or ppis)
	bool ppippi_ppippi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool pxpx_pypy_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool ppippi_pzpz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool pzpz_ppippi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool pzpz_pzpz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	bool spz_ppippi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (spx or pxs)
	bool ppippi_spz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (spx or pxs)
	bool spz_pzpz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (spx or pxs)
	bool pzpz_spz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (spx or pxs)
	bool spz_spz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (spx or pxs)
	bool sppi_ppipz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (sppi or ppis) and (ppipx or pxppi)
	bool ppipz_sppi_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (sppi or ppis) and (ppipx or pxppi)
	bool ppipz_ppipz_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (ppipx or pxpi)
	bool pxpy_pxpy_condition(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma); // (pypz or pzpy)


public:
	int check_two_center_type(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);
	int check_one_center_type(int orbital_mi, int orbital_ni, int orbital_lambda, int orbital_sigma);

};

#endif
