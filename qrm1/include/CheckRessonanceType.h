#ifndef CHECKRESSONANCETYPE_H
#define CHECKRESSONANCETYPE_H

class CheckRessonanceType
{
private:
	bool check_1s_1s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_1s_2s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2s_1s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_1s_2pz(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2pz_1s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2s_2s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2s_2pz(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2pz_2s(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2pz_2pz(int q_number_A, int q_number_B, int mi, int ni);
	bool check_2ppi_2ppi(int q_number_A, int q_number_B, int mi, int ni);
	bool check_s_pi(int mi, int ni); //always zero
	bool check_p1_p2(int mi, int ni); //always zero


public:
	int choose_ressonance_expression(int q_number_A, int q_number_B, int mi, int ni);


};

#endif
