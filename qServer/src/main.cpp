#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "Socket.h"
#include "ClientSocket.h"
#include "SocketException.h"

using namespace std;

vector<bool> done;
vector<int> waysToRun;

bool is_file_exist(const char *fileName)
{
	ifstream infile(fileName);
	return infile.good();
}

string getSystem(int nSystem)
{
	switch (nSystem)
	{
	case 0:
		return "qinputh2p.txt";
	case 1:
		return "qinputh2.txt";
	case 2:
		return "qinputh3p.txt";
	case 3:
		return "qinputh4p.txt";
	case 4:
		return "qinputh5p.txt";
	case 5:
		return "qinputh6p.txt";
	case 6:
		return "qinputh7p.txt";
	case 7:
		return "qinputh8p.txt";
	case 8:
		return "qinputh9p.txt";
	default:
		cout << "error on - string setSystem(int nSystem)" << endl
			<< "contact developers" << endl;
		exit(7);
	}
}

int findSpotToRun(int iRun)
{
	int size = waysToRun.size();
	for (int i = 0; i<size; i++)
	{
		if (waysToRun[i] == -1)
		{
			waysToRun[i] = iRun;
			return i;
		}
	}

	cout << "Erro em findSpotToRun " << endl
		<< "Isso nao devia acontecer, contate desenvolvedores"
		<< endl;
	return -1;
}

void freeWayToRun(int iRun)
{
	int size = waysToRun.size();
	for (int i = 0; i<size; i++)
	{
		if (waysToRun[i] == iRun)
		{
			waysToRun[i] = -1;
			return;
		}
	}
	//cout << "Erro em freeWayToRun " << endl
	//	<< "Isso nao devia acontecer, contate desenvolvedores"
	//	<< endl;
}
// nas primeiras ele tenta chegar nesse ponto antes mesmo de ter comecado a rodar.
// o ssh e lento.

// ssh fred@150.161.5.110 /home/fred/rindo/rindo/client/teste2/rindo_client.exe 46
void chooseRunMethod(int way, int i)
{
	string run, pega;
	switch (way)
	{
	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
		run = "./qRM1.exe " + getSystem(i);
		system(run.c_str());
		break;

	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
		//		run = clientPath2 + "rindo_client.exe  " + to_string(i);
		//		system(run.c_str());
		//		pega = "sftp -b ./batch/" + inpNames[i] + ".b " + clientName2;
		//		system(pega.c_str());
		//nova forma
		break;
	case 14:
	case 15:
	case 16:
	case 17:
	case 18:
	case 19:
		//		run = clientPath1 + "rindo_client.exe  " + to_string(i);
		//		system(run.c_str());
		break;
	default:
		cout << "Erro em runRindo" << endl
			<< "Isso nao devia acontecer, contate desenvolvedores"
			<< endl;
		exit(1);
	}
}

void rodaProcesso(int i)
{
	chooseRunMethod(findSpotToRun(i), i);
	//cout << "comecou o :  " << i << endl;
	done[i] = true;
	//cout << i << " finalizado" << endl;
}


int main()
{
	string fitness;
	int repete = 0;
	if (!is_file_exist("restartEnd.ga"))
	{
		while(repete==0)
		{
			try
			{
				ClientSocket client_socket("localhost", 30000);
				try
				{
					client_socket << "r";
					client_socket >> fitness;
					cout << "pegou:   " << fitness << endl;
					repete=1;
				}
				catch (SocketException&) {}
			}
			catch (SocketException& e)
			{
				std::cout << "Exception was caught:" << e.description() << "\n";
				repete = 0;
				sleep(0.001);
			}
		}
	}
	else
	{


		const int totalJobs = 5;

		ifstream iproc_("proc.txt");
		int nproc;
		iproc_ >> nproc;
		iproc_.close();

		done.resize(totalJobs);
		vector<bool> activeProcess(totalJobs);
		for (int i = 0; i < totalJobs; i++)
		{
			done[i] = false;
			activeProcess[i] = false;
		}

		thread t[totalJobs];
		waysToRun.resize(nproc);
		for (int i = 0; i < nproc; i++)
		{
			waysToRun[i] = -1;
			t[i] = thread(rodaProcesso, i);
			activeProcess[i] = true;
		}

		//loop pra enfiar processos onde for necessario
		bool finished = false;
		while (!finished)
		{
			for (int i = 0; i < totalJobs; i++)
			{
				//cout << "done:  " << done[i] << "  active " << activeProcess[i] << endl;
				//cin.get();

				if (done[i] && activeProcess[i])
				{
					t[i].join();
					freeWayToRun(i);
					activeProcess[i] = false;

					if (i == (totalJobs - 1))
					{
						//cout << "acabou" << endl;
						finished = true;
						break;
					}
					else
					{
						for (int ii = i + 1; ii < totalJobs; ii++)
						{
							if ((!activeProcess[ii]) && (!done[ii]))
							{
								t[ii] = thread(rodaProcesso, ii);
								activeProcess[ii] = true;
								break;
							}
						}
					}
				}
			}
		}


		for (int i = 0; i < totalJobs; i++)
		{
			if (t[i].joinable())
			{
				t[i].join();
			}
		}

		ifstream in_;
		double error = 0.0e0;
		double aux;
		in_.open("h2p.ga");
		in_ >> aux;
		error += aux;
		in_.close();
		in_.open("h2.ga");
		in_ >> aux;
		error += aux;
		in_.close();
		in_.open("h3p.ga");
		in_ >> aux;
		error += aux;
		in_.close();
		in_.open("h4p.ga");
		in_ >> aux;
		error += aux;
		in_.close();
		in_.open("h5p.ga");
		in_ >> aux;
		error += aux;
		in_.close();

		stringstream convert;
		convert << error / 5.0e0;
		convert >> fitness;

		ofstream restart_;
		restart_.open("restart.ga", ios_base::app);
		restart_ << error / 5.0e0 << endl;
		restart_.close();
	}

	ofstream of_("fitness.ga");
	of_  << fitness << endl;
	of_.close();

	return 0;
}

/* teste pra ver se e deterministico
	ofstream gen_;
	gen_.open("genes.ga", ios_base::app);
	double auxHpar;
	ifstream hpar_("hparam.txt");
	for (int ii = 0; ii < 13; ii++)
	{
		hpar_ >> auxHpar;
		gen_ << auxHpar << "   ";
			
	}
	hpar_.close();
	gen_ << endl;
	gen_.close();
*/

