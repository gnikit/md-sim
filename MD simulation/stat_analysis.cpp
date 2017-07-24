#include "stat_analysis.h"

Stat_Analysis::Stat_Analysis(std::string PATH, vec1d A_LIST) {
	// Do everything in the constructor
	// No need for a Run function
	_path = PATH;
	_A_list = A_LIST;
}

Stat_Analysis::~Stat_Analysis() {
}

void Stat_Analysis::ReadFromFile(std::vector<double> &x, const std::string &file_name) {
	/*
	Reads from a stream that already exists for a file that is already placed in the
	directory and appends the data into a 1D vector.
	*/
	std::ifstream read_file(file_name);
	assert(read_file.is_open());

	std::copy(std::istream_iterator<double>(read_file),
		std::istream_iterator<double>(), std::back_inserter(x));

	read_file.close();
}

double Stat_Analysis::Mean(std::vector<double> &x) {
	size_t n = 0;
	double temp = 0;
	for (size_t i = 0; i < x.size(); i++) {
		temp += x.at(i);
		n++;
	}
	temp /= n; // this is the mean
	m1 = temp;

	return m1;
}

void Stat_Analysis::StaticDataProcessing(size_t n)
/*
Used to average the energies and pressures for all
values of A for a specific number n of the potential
strength
*/
{
	std::vector<double> A{ 0, 0.25, 0.5,  0.75, 1, 1.25, 1.5, 1.75,
		2, 2.25, 2.50, 2.75, 3, 3.5,  4 };
	std::vector<double> temp;

	std::ofstream data;
	std::string K, U, Tot, Pc, path, sep, power, a, txt;
	path = "../../Archives of Data/Density 0.5/Isothermal~step 5000/"; // warning change step count
	sep = "~";														// TODO: generic step change
	txt = ".txt";
	// path = "\0";
	power = std::to_string(n);
	std::string name = path + "data" + power + txt;
	data.open(name, std::ios::out | std::ios::trunc);

	// data << "Power:\t" << n << endl;
	data << "# A:\tK:\tU:\tTot:\tPc:" << endl;

	for (size_t i = 0; i < A.size(); i++) {
		K.clear();
		U.clear();
		Tot.clear();
		Pc.clear();
		temp1.clear();
		temp2.clear();
		temp3.clear();
		temp4.clear();
		K = "KinEn";
		U = "PotEn";
		Tot = "TotEn";
		Pc = "PressureC";
		std::stringstream stream;
		stream << std::fixed << std::setprecision(2) << A.at(i);
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