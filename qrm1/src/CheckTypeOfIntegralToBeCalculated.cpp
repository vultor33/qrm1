#include "CheckTypeOfIntegralToBeCalculated.h"

#include <iostream>
#include <stdlib.h>

using namespace std;

int CheckTypeOfIntegralToBeCalculated::check_one_center_type(int mi, int ni, int lambda, int sigma)
{
	if (ss_ss_condition(mi, ni, lambda, sigma))// (ss|ss)
	{
		return 1;
	}
	else if (ss_pp_condition(mi, ni, lambda, sigma))
	{
		return 2;
	}
	else if (sp_sp_condition(mi, ni, lambda, sigma))
	{
		return 3;
	}
	else if (pp_pp_condition(mi, ni, lambda, sigma))
	{
		return 4;
	}
	else if (pp_plpl_condition(mi, ni, lambda, sigma))
	{
		return 5;
	}
	else if (ppl_ppl_condition(mi, ni, lambda, sigma))
	{
		return 6;
	}
	else
	{
		return 0;// atencao!!!!!
		cout << "integral nao cadastrada - consultar rotina check_one_center_type" << endl;
		exit(4);
	}
	return 0;
}

bool CheckTypeOfIntegralToBeCalculated::is_p_orbital(int orb)
{
	//true if orb == p
	return ((orb == 1) || (orb == 2) || (orb == 3));
}

bool CheckTypeOfIntegralToBeCalculated::ss_ss_condition(int mi, int ni, int lambda, int sigma)
{
	return (mi == 0) && (ni == 0) && (lambda == 0) && (sigma == 0);
}

bool CheckTypeOfIntegralToBeCalculated::sp_sp_condition(int mi, int ni, int lambda, int sigma)
{

	return 
		((mi == 0) && (is_p_orbital(ni)) && (lambda == 0) && (is_p_orbital(sigma))) //spsp
		||
		((mi == 0) && (is_p_orbital(ni)) && (sigma == 0) && (is_p_orbital(lambda))) //spps
		||
		((ni == 0) && (is_p_orbital(mi)) && (lambda == 0) && (is_p_orbital(sigma))) //pssp
		||
		((ni == 0) && (is_p_orbital(mi)) && (sigma == 0) && (is_p_orbital(lambda))) //psps
		;
}

bool CheckTypeOfIntegralToBeCalculated::ss_pp_condition(int mi, int ni, int lambda, int sigma)
{
	return 
		((mi == 0) && (ni == 0) && is_p_orbital(lambda) && (lambda==sigma)) //sspp
		||
		((lambda == 0) && (sigma == 0) && is_p_orbital(mi) && (mi == ni)) //ppss
		;
}

bool CheckTypeOfIntegralToBeCalculated::pp_pp_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == ni) && (ni == sigma)&&(sigma == lambda)) //all equal
		&&
		is_p_orbital(mi);  		
}

bool CheckTypeOfIntegralToBeCalculated::pp_plpl_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == ni) && (is_p_orbital(mi)))//pp
		&&
		(mi != lambda)
		&&
		((lambda == sigma) && (is_p_orbital(lambda)))//plpl
		);
}

bool CheckTypeOfIntegralToBeCalculated::ppl_ppl_condition(int mi, int ni, int lambda, int sigma)
{

	return
		((mi!=ni) && (mi == lambda) && (is_p_orbital(mi)) && (ni == sigma) && (is_p_orbital(ni))) //ppl_ppl and plp_plp
		||
		((mi != ni) && (mi == sigma) && (is_p_orbital(mi)) && (ni == lambda) && (is_p_orbital(ni))) //ppl_plp and plp_ppl
		;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TWO CENTER ROUTINES //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Aqui vao sair numeros que vao escolher as expressoes que serao calculadas
int CheckTypeOfIntegralToBeCalculated::check_two_center_type(int mi, int ni, int lambda, int sigma)
{
	//cout << "mi: " << mi << "  ni  " << ni << " lam  " << lambda << "  sig  " << sigma << endl;
	if (ss_ss_condition(mi, ni, lambda, sigma))// (ss|ss)
	{
		return 1;
	}
	else if (ss_spz_condition(mi, ni, lambda, sigma))
	{
		return 2;
	}
	else if (spz_ss_condition(mi, ni, lambda, sigma))
	{
		return 3;
	}
	else if (ss_pzpz_condition(mi, ni, lambda, sigma))
	{
		return 4;
	}
	else if (pzpz_ss_condition(mi, ni, lambda, sigma))
	{
		return 5;
	}
	else if (ss_ppippi_condition(mi, ni, lambda, sigma))
	{
		return 6;
	}
	else if (ppippi_ss_condition(mi, ni, lambda, sigma))
	{
		return 7;
	}
	else if (sppi_sppi_condition(mi, ni, lambda, sigma))
	{
		return 8;
	}
	else if (ppippi_ppippi_condition(mi, ni, lambda, sigma))
	{
		return 9;
	}
	else if (pxpx_pypy_condition(mi, ni, lambda, sigma))
	{
		return 10;
	}
	else if (ppippi_pzpz_condition(mi, ni, lambda, sigma))
	{
		return 11;
	}
	else if (pzpz_ppippi_condition(mi, ni, lambda, sigma))
	{
		return 12;
	}
	else if (pzpz_pzpz_condition(mi, ni, lambda, sigma))
	{
		return 13;
	}
	else if (spz_ppippi_condition(mi, ni, lambda, sigma))
	{
		return 14;
	}
	else if (ppippi_spz_condition(mi, ni, lambda, sigma))
	{
		return 15;
	}
	else if (spz_pzpz_condition(mi, ni, lambda, sigma))
	{
		return 16;
	}
	else if (pzpz_spz_condition(mi, ni, lambda, sigma))
	{
		return 17;
	}
	else if (spz_spz_condition(mi, ni, lambda, sigma))
	{
		return 18;
	}
	else if (sppi_ppipz_condition(mi, ni, lambda, sigma))
	{
		return 19;
	}
	else if (ppipz_sppi_condition(mi, ni, lambda, sigma))
	{
		return 20;
	}
	else if (ppipz_ppipz_condition(mi, ni, lambda, sigma))
	{
		return 21;
	}
	else if (pxpy_pxpy_condition(mi, ni, lambda, sigma))
	{
		return 22;
	}
	else
	{
		return 99;
		cout << "integral nao cadastrada - consultar rotina check_integral_type" << endl;
		exit(4);
	}
	return 0;
}

bool CheckTypeOfIntegralToBeCalculated::ss_spz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 0) && (ni == 0)) //ss
		&&
		(
		(lambda == 0) && (sigma == 3)
		||
		(lambda == 3) && (sigma == 0)
		); //spx or pxs		
}

bool CheckTypeOfIntegralToBeCalculated::spz_ss_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		(mi == 0) && (ni == 3)
		||
		(mi == 3) && (ni == 0)
		)//spx or pxs
		&&
		((lambda == 0) && (sigma == 0)) //ss
		; 		
}

bool CheckTypeOfIntegralToBeCalculated::ss_pzpz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 0) && (ni == 0)) //ss
		&&
		((lambda == 3) && (sigma == 3)); //pxpx 		
}

bool CheckTypeOfIntegralToBeCalculated::pzpz_ss_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 3) && (ni == 3)) //pxpx
		&&
		((lambda == 0) && (sigma == 0)); //ss
}

bool CheckTypeOfIntegralToBeCalculated::ss_ppippi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 0) && (ni == 0)) //ss
		&&
		((lambda == sigma) && ((sigma == 2)||(sigma==1))); //pxpx 		
}

bool CheckTypeOfIntegralToBeCalculated::ppippi_ss_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((lambda == 0) && (sigma == 0)) //ss
		&&
		((mi == ni) && ((mi == 2) || (mi == 1))); //pxpx 		
}


bool CheckTypeOfIntegralToBeCalculated::sppi_sppi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 0) && ((ni == 2) || (ni == 1)) && (ni == sigma) && (lambda==0))
		||
		((ni == 0) && ((mi == 2) || (mi == 1)) && (mi == sigma) && (lambda==0))
		||
		((mi == 0) && ((ni == 2) || (ni == 1)) && (ni == lambda) && (sigma==0))
		||
		((ni == 0) && ((mi == 2) || (mi == 1)) && (mi == lambda) && (sigma==0))
		;
}

bool CheckTypeOfIntegralToBeCalculated::ppippi_ppippi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == ni) && (mi == lambda) && (mi==sigma))&&
		((mi==2)||(mi==1))
		);
}

bool CheckTypeOfIntegralToBeCalculated::pxpx_pypy_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == ni) && (mi == 2)) &&
		((lambda == sigma) && (lambda == 1))
		)||
		(
		((mi == ni) && (mi == 1)) &&
		((lambda == sigma) && (lambda == 2))
		);
}

bool CheckTypeOfIntegralToBeCalculated::ppippi_pzpz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == ni) && ((mi == 2) || (mi==1)))&&
		((lambda == sigma) && (lambda == 3));
}

bool CheckTypeOfIntegralToBeCalculated::pzpz_ppippi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((lambda == sigma) && ((lambda == 2) || (sigma == 1))) &&
		((mi == ni) && (mi == 3));
}

bool CheckTypeOfIntegralToBeCalculated::pzpz_pzpz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == ni) && (mi == lambda) && (mi == sigma) && (mi == 3));
}

bool CheckTypeOfIntegralToBeCalculated::spz_ppippi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == 0) && (ni == 3)) || ((mi == 3) && (ni == 0))
		) &&
		((lambda == sigma) && ((lambda == 2) || (lambda == 1)));
}

bool CheckTypeOfIntegralToBeCalculated::ppippi_spz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((lambda == 0) && (sigma == 3)) || ((lambda == 3) && (sigma == 0))
		) &&
		((mi == ni) && ((mi == 2) || (mi == 1)));
}

bool CheckTypeOfIntegralToBeCalculated::spz_pzpz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == 0) && (ni == 3)) || ((mi == 3) && (ni == 0))
		) &&
		((lambda == sigma) && (lambda == 3));
}

bool CheckTypeOfIntegralToBeCalculated::pzpz_spz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((lambda == 0) && (sigma == 3)) || ((lambda == 3) && (sigma == 0))
		) &&
		((mi == ni) && (mi == 3));
}

bool CheckTypeOfIntegralToBeCalculated::spz_spz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		(
		((mi == 0) && (ni == 3)) || ((mi == 3) && (ni == 0))
		) &&
		((lambda == 0) && (sigma == 3)) || ((lambda == 3) && (sigma == 0));
}

bool CheckTypeOfIntegralToBeCalculated::sppi_ppipz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 0) && ((ni == 2)||(ni==1)) && (lambda == ni) && (sigma == 3))//(sppi|ppipx)
		||
		((ni == 0) && ((mi == 2)||(mi==1)) && (lambda == mi) && (sigma == 3))//(ppis|ppipx)
		||
		((mi == 0) && ((ni == 2) || (ni == 1)) && (sigma == ni) && (lambda == 3))//(sppi|pxppi)
		||
		((ni == 0) && ((mi == 2) || (mi == 1)) && (sigma == mi) && (lambda == 3))//(ppis|pxppi)
		;
}

bool CheckTypeOfIntegralToBeCalculated::ppipz_sppi_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((lambda == 0) && ((sigma == 2) || (sigma == 1)) && (mi == sigma) && (ni == 3))//(ppipx|sppi)
		||
		((sigma == 0) && ((lambda == 2) || (lambda == 1)) && (mi == lambda) && (ni == 3))//(ppipx|ppis)
		||
		((lambda == 0) && ((sigma == 2) || (sigma == 1)) && (ni == sigma) && (mi == 3))//(pxppi|sppi)
		||
		((sigma == 0) && ((lambda == 2) || (lambda == 1)) && (ni == lambda) && (mi == 3))//(pxppi|ppi)
		;
}

bool CheckTypeOfIntegralToBeCalculated::ppipz_ppipz_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 3) && ((ni == 2) || (ni == 1)) && (lambda == ni) && (sigma == 3))//(pxppi|ppipx)
		||
		((ni == 3) && ((mi == 2) || (mi == 1)) && (lambda == mi) && (sigma == 3))//(ppipx|ppipx)
		||
		((mi == 3) && ((ni == 2) || (ni == 1)) && (sigma == ni) && (lambda == 3))//(pxppi|pxppi)
		||
		((ni == 3) && ((mi == 2) || (mi == 1)) && (sigma == mi) && (lambda == 3))//(ppipx|pxppi)
		;
}

bool CheckTypeOfIntegralToBeCalculated::pxpy_pxpy_condition(int mi, int ni, int lambda, int sigma)
{
	return
		((mi == 2) && (ni == 1) && (lambda == 2) && (sigma == 1))
		||
		((ni == 2) && (mi == 1) && (lambda == 2) && (sigma == 1))
		||
		((mi == 2) && (ni == 1) && (sigma == 2) && (lambda == 1))
		||
		((ni == 2) && (mi == 1) && (sigma == 2) && (lambda == 1))
		;
}
