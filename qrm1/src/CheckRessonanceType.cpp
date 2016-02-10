#include "CheckRessonanceType.h"

#include <iostream>
#include <stdlib.h>

using namespace std;


int CheckRessonanceType::choose_ressonance_expression(int q_number_A, int q_number_B, int mi, int ni)
{
	if (check_1s_1s(q_number_A,q_number_B,mi,ni))
	{
		return 1;
	}
	else if (check_1s_2s(q_number_A, q_number_B, mi, ni))
	{
		return 2;
	}
	else if (check_2s_1s(q_number_A, q_number_B, mi, ni))
	{
		return 3;
	}
	else if (check_1s_2pz(q_number_A, q_number_B, mi, ni))
	{
		return 4;
	}
	else if (check_2pz_1s(q_number_A, q_number_B, mi, ni))
	{
		return 5;
	}
	else if (check_s_pi(mi, ni))
	{
		return 0;
	}
	else if (check_2s_2s(q_number_A, q_number_B, mi, ni))
	{
		return 6;
	}
	else if (check_2s_2pz(q_number_A, q_number_B, mi, ni))
	{
		return 7;
	}
	else if (check_2pz_2s(q_number_A, q_number_B, mi, ni))
	{
		return 8;
	}
	else if (check_2pz_2pz(q_number_A, q_number_B, mi, ni))
	{
		return 9;
	}
	else if (check_2ppi_2ppi(q_number_A, q_number_B, mi, ni))
	{
		return 10;
	}
	else if (check_p1_p2(mi, ni))
	{
		return 0;
	}
	else
	{
		cout << " checando tipos ressonancia: ------ nao cadastrado" << endl;
		return 99;
		cout << "base nao cadastrada - erro em choose_ressonance_expression" << endl;
		exit(4);
	}
}

bool CheckRessonanceType::check_1s_1s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 1) &&
		(q_number_B == 1) &&
		(mi == 0) &&
		(ni == 0);
}

bool CheckRessonanceType::check_1s_2s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 1) &&
		(q_number_B == 2) &&
		(mi == 0) &&
		(ni == 0);
}

bool CheckRessonanceType::check_2s_1s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 1) &&
		(mi == 0) &&
		(ni == 0);
}

bool CheckRessonanceType::check_1s_2pz(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 1) &&
		(q_number_B == 2) &&
		(mi == 0) &&
		(ni == 3);
}

bool CheckRessonanceType::check_2pz_1s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 1) &&
		(mi == 3) &&
		(ni == 0);
}

bool CheckRessonanceType::check_2s_2s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 2) &&
		(mi == 0) &&
		(ni == 0);
}

bool CheckRessonanceType::check_s_pi(int mi, int ni)
{
	return
		((mi == 0) && ((ni == 2)||(ni == 1)))
		||
		((ni == 0) && ((mi == 2) || (mi == 1)))
		;
}

bool CheckRessonanceType::check_2s_2pz(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 2) &&
		(mi == 0) &&
		(ni == 3);
}

bool CheckRessonanceType::check_2pz_2s(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 2) &&
		(mi == 3) &&
		(ni == 0);
}

bool CheckRessonanceType::check_2pz_2pz(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 2) &&
		(mi == 3) &&
		(ni == 3);
}

bool CheckRessonanceType::check_2ppi_2ppi(int q_number_A, int q_number_B, int mi, int ni)
{
	return
		(q_number_A == 2) &&
		(q_number_B == 2) &&
		(mi == ni) &&
		((ni == 2)||(ni==1));
}


bool CheckRessonanceType::check_p1_p2(int mi, int ni)
{
	return
		(mi != ni) &&
		((mi==1) || (mi == 2) || (mi == 3))&&
		((ni == 1) || (ni == 2) || (ni == 3))
		;
}







