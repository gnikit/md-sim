#include "MD.h"
#include <cstdio>

using namespace std;

FILE *stream; 
vector<double> temp1, temp2, temp3, temp4;
double m1;



void ReadFromFile(std::vector<double> &x, const std::string &file_name);
double Mean(std::vector<double> &x);
void StaticDataProcessing(size_t n);


int main()
{
	//freopen_s(&stream, "LOG.txt", "a+", stderr);
	size_t steps = 5001;
	vector<double> A_parameter = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25, 2.50, 2.75, 3.0, 3.5, 4.0 };
	vector<size_t> power = { 6,7,8,9,10,12 };
	//cerr << "Log file of the MD simulation" << endl;
	//cerr << "Power:\tA:\tsteps:" << endl;
	size_t num = 1;
	char full = 'y';
	//cout << "Run full simulation?"; cin >> full;
	// Currently this code is not executing
	//if (full == 'y')
	//{
	//	for (size_t n = 0; n < power.size(); n++)
	//	{
	//		for (size_t a = 0; a < A_parameter.size(); a++)
	//		{
	//			cout << "run: " << num << endl;
	//			MD run(power.at(n), A_parameter.at(a), steps);
	//			run.Simulation();
	//			//cerr << power.at(n) << "\t" << A_parameter.at(a) << "\t" << steps << "rho = 0.8" << endl;
	//			num++;
	//		}
	//	}
	//}

	//vector<size_t> step = { 15000 };
	//
	//if (full == 'y')
	//{
	//	for (size_t n = 0; n < step.size(); n++)
	//	{
	//		for (size_t a = 0; a < A_parameter.size(); a++)
	//		{
	//			cout << "run: " << num << endl;
	//			MD run(6, A_parameter.at(a), step.at(n));
	//			run.Simulation();
	//			cerr << 6 << '\t' << A_parameter.at(a) << '\t' << step.at(n) << endl;
	//			num++;
	//		}
	//	}
	//}
	//StaticDataProcessing(6);


	//for (size_t j = 0; j < power.size(); j++)
	//{
	//	MD run(power.at(j), 1.1, 5001);
	//	run.Simulation();
	//}
	//StaticDataProcessing(11);

	//MD run(20, 1.0001, steps);
	//run.Simulation();
	//MD run1(11, 1.1, steps);
	//run1.Simulation();
	//
	//MD run1(6, 0.5, steps);
	//run1.Simulation(); cerr << "U: 6\tA:0.5" << endl;
	//
	//MD run2(6, 0.75, steps);
	//run2.Simulation(); cerr << "U: 6\tA:0.75" << endl;
	//
	//MD run3(6, 1, steps);
	//run3.Simulation(); cerr << "U: 6\tA:1" << endl;
	//
	//MD run4(6, 1.25, steps);
	//run4.Simulation(); cerr << "U: 6\tA:1.25" << endl;
	//
	//MD run5(6, 1.5, steps);
	//run5.Simulation(); cerr << "U: 6\tA:1.5" << endl;
	
	// Static Data Processing
	//for (size_t i = 6; i < 13; i++)
	//{	
	//	if (i == 11)
	//	{	i++;	}
	//	StaticDataProcessing(i);
	//}

	system("pause");

}

void ReadFromFile(std::vector<double> &x, const std::string &file_name)
/*
Reads from a stream that already exists for a file that is already placed in the directory
and appends the data into a 1D vector.
*/
{
	std::ifstream read_file(file_name);
	assert(read_file.is_open());

	std::copy(std::istream_iterator<double>(read_file), std::istream_iterator<double>(),
		std::back_inserter(x));

	read_file.close();
}

double Mean(std::vector<double> &x)
{
	size_t n = 0;
	double temp = 0;
	for (size_t i = 0; i < x.size(); i++)
	{
		temp += x.at(i);
		n++;
	}
	temp /= n;	// this is the mean
	m1 = temp;

	return m1;
}

void StaticDataProcessing(size_t n)
{
	using namespace std; 

	vector<double> A{ 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.50, 2.75, 3, 3.5, 4 };
	vector<double> temp;

	ofstream data;
	string K, U, Tot, Pc, path, sep, power, a, txt;
	path = "../../Archives of Data/Density 0.5/Isothermal~step 5000/"; sep = "~"; txt = ".txt";
	//path = "\0";
	power = std::to_string(n);
	string name = path + "data" + power + txt;
	data.open(name, std::ios::out | std::ios::trunc);

	//data << "Power:\t" << n << endl;
	data << "#A:\tK:\tU:\tTot:\tPc:" << endl;

	for (size_t i = 0; i < A.size(); i++)
	{
		K.clear(); U.clear();  Tot.clear(); Pc.clear();
		temp1.clear(); temp2.clear(); temp3.clear(); temp4.clear();
		K = "KinEn"; U = "PotEn"; Tot = "TotEn"; Pc = "PressureC";
		std::stringstream stream;
		stream << fixed << setprecision(2) << A.at(i);
		a = stream.str();
		K = path + K + power + sep + a + txt;
		U = path + U + power + sep + a + txt;
		Tot = path + Tot + power + sep + a + txt;
		Pc = path + Pc + power + sep + a + txt;

		ReadFromFile(temp1, K); 
		ReadFromFile(temp2, U);
		ReadFromFile(temp3, Tot);
		ReadFromFile(temp4, Pc);

		Mean(temp1);
		data.precision(5);
		data << A.at(i) << '\t' << m1 << '\t';
		Mean(temp2); 
		data << m1 << '\t';
		Mean(temp3);
		data << m1 << '\t';
		Mean(temp4);
		data << m1 << endl;
	}
}
