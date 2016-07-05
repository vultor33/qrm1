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
int nproc;
int nhydrogens;

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
}


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
	case 14:
	case 15:
	case 16:
	case 17:
	case 18:
	case 19:
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



bool qmodel()
{
// set model here without input
// double precision
	remove("hparam.txt");
	ifstream inputServer_("inputServer.txt");
	int restart, model;
	inputServer_ >> restart;
	inputServer_ >> model;
	inputServer_ >> nproc;
	inputServer_ >> nhydrogens;
	inputServer_.close();


	ifstream points_("point.txt");
	int nPoints;
	vector<double> points; 
	points_ >> nPoints;
	points.resize(nPoints);
	for(int i=0; i<nPoints; i++)
	{
		points_ >> points[i];
	}
	points_.close();

	// models:  0 - qrm1 ; 1 - qover ; 2 - qint; 3 - qoverqint ; 4 - qalfa; 5 - qoverqalfa
	// 1-over; 2-int; 3-alfa; 4-gauss
	ofstream hparam_("hparam.txt");
	if((model == 8)||(model==9))
		hparam_ << 0 << endl;
	else
		hparam_ << 5 << endl;
	
	if(model < 6)
	{
		for(int i=0; i<8; i++)
			hparam_ << setprecision(16) << points[i] << endl;
	}
	switch(model)
	{
		case 0:
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 1:
			hparam_ << setprecision(16) << points[8] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 2:
			hparam_ << 1.0e0 << endl;
			hparam_ << setprecision(16) << points[8] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 3:
			hparam_ << setprecision(16) << points[8] << endl;
			hparam_ << setprecision(16) << points[9] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 4:
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << setprecision(16) << points[8] << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 5:
			hparam_ << setprecision(16) << points[8] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << setprecision(16) << points[9] << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 6:
			hparam_ << 1.1901  << endl;
			hparam_ << 1.2784  << endl;
			hparam_ << 0.9300 << endl;
			hparam_ << 1.4933 << endl;
			hparam_ << setprecision(16) << points[0] << endl;
			hparam_ << 1.1373 << endl;
			hparam_ << 2.4347 << endl;
			hparam_ << 1.1293 << endl;
			hparam_ << setprecision(16) << points[1] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 7:
			hparam_ << 1.1972  << endl;
			hparam_ << 1.3021  << endl;
			hparam_ << 0.9381 << endl;
			hparam_ << 1.4816 << endl;
			hparam_ << setprecision(16) << points[0] << endl;
			hparam_ << 0.6860 << endl;
			hparam_ << 12.1854 << endl;
			hparam_ << 1.7997 << endl;
			hparam_ << setprecision(16) << points[1] << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;
		
		case 8:
			hparam_ << 1.1901  << endl;
			hparam_ << 1.2784  << endl;
			hparam_ << 0.9300 << endl;
			hparam_ << 1.4933 << endl;
			hparam_ << setprecision(16) << points[0] << endl;
			hparam_ << 1.1373 << endl;
			hparam_ << 2.4347 << endl;
			hparam_ << 1.1293 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		case 9:
			hparam_ << 1.1972  << endl;
			hparam_ << 1.3021  << endl;
			hparam_ << 0.9381 << endl;
			hparam_ << 1.4816 << endl;
			hparam_ << setprecision(16) << points[0] << endl;
			hparam_ << 0.6860 << endl;
			hparam_ << 12.1854 << endl;
			hparam_ << 1.7997 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			hparam_ << 1.0e0 << endl;
			break;

		default:
			cout << "ERRO EM qServer.x:  qmodel()" << endl;
			exit(1);
	}

	if(restart==0)
		return true;
	else
		return false;
}

int main()
{
	string fitness;
	int repete = 0;
	if (qmodel())
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
		int totalJobs = nhydrogens;

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
		switch(nhydrogens)
		{
			case 9:
				in_.open("h9p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 8:
				in_.open("h8p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 7:
				in_.open("h7p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 6:
				in_.open("h6p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 5:
				in_.open("h5p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 4:
				in_.open("h4p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 3:
				in_.open("h3p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 2:
				in_.open("h2.ga");
				in_ >> aux;
				error += aux;
				in_.close();
			case 1:
				in_.open("h2p.ga");
				in_ >> aux;
				error += aux;
				in_.close();
				break;
			default:
				cout << "ERROR ON getFitness.x" << endl;
				exit(1);
		}	
		
		stringstream convert;
		convert << setprecision(16) << error / nhydrogens;
		convert >> fitness;

		//ofstream restart_;
		//restart_.open("restart.fit", ios_base::app);
		//restart_ << fitness << endl;
		//restart_.close();
	}

	ofstream of_("fitness.txt");
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

