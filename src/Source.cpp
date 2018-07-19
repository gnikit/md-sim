#include "MD.h"
#include "isomorph.h"
#include "stat_analysis.h"
#include <thread>
#include <string>


#define STEPS 100	//10000
#define PARTICLES 1000 //1000
typedef std::vector<double> vec1d;

/* Windows working directory for Archives of Data */
std::string dir_windows = "C:/Code/C++/MD-simulation/Archives of Data/testing/gaussian/long_quench/";
/* Linux working directory */
std::string dir_linux = "/home/gn/Desktop/test_data/delete/";

void MakeDataBase();

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
	
	vec1d temp = LinearSpacedArray(0.001, 0.01, 10);
	for (size_t i = 0; i < temp.size(); i++) {
		MD* run1 = new MD(dir_linux, STEPS, false);		// fluid compression set to false
		MD* run2 = new MD(dir_linux, STEPS, true);
		MD* run3 = new MD(dir_linux, STEPS, true);
		MD* run4 = new MD(dir_linux, STEPS, true);

		std::thread th1(&MD::Simulation, run1, 0.05, 0.5, 12, 0.0);
		/* std::thread th2(&MD::Simulation, run2, 0.05, 0.5, 12, 0.5);
		std::thread th3(&MD::Simulation, run3, 0.05, 0.5, 12, 0.1);
		std::thread th4(&MD::Simulation, run4, 0.05, 0.5, 12, 1.5);
 */
		th1.join();/* th2.join(); th3.join(); th4.join(); */
		delete run1, run2, run3, run4;
	}
}

void MakeDataBase() {
	/*
	* Builds a relatively complete database mapping the rho,T,n,a space
	* of te MD fluid. Performs Also statistical averaging at the end with
	* quantities such as Pc, U, K, E, Pk.
	* NOTE: It takes a long time to complete and generates multiple files!
	*/
	size_t num = 1;
	std::vector<size_t> n = { 6, 8, 10, 12 };
	std::vector<double> rho = { 0.5, 1.0, 1.5, 2.0 };
	std::vector<double> T = { 0.5, 1.0, 1.5, 2.0 };
	std::vector<double> A1 = LinearSpacedArray(0, 1, 5);
	std::vector<double> A2 = LinearSpacedArray(1.25, 2.25, 5);
	std::vector<double> A3 = LinearSpacedArray(2.50, 4.50, 5);
	std::vector<double> A4 = LinearSpacedArray(0.6, 1.0, 5);  // Includes a=1.0 again    

	for (size_t d = 0; d < rho.size(); d++) {
		for (size_t t = 0; t < T.size(); t++) {
			for (size_t i = 0; i < n.size(); i++) {
				for (size_t j = 0; j < A1.size(); j++) {
					std::cout << "rho: " << rho[d] << " T: " << T[t] <<
						"n: " << n[i] << " A: " << A1[j] << " run num: " << num << std::endl;

					MD* run1 = new MD(dir_windows, STEPS);
					MD* run2 = new MD(dir_windows, STEPS);
					MD* run3 = new MD(dir_windows, STEPS);
					MD* run4 = new MD(dir_windows, STEPS);

					std::thread th1(&MD::Simulation, run1, rho[d], T[t], n[i], A1[j]);
					std::thread th2(&MD::Simulation, run2, rho[d], T[t], n[i], A2[j]);
					std::thread th3(&MD::Simulation, run3, rho[d], T[t], n[i], A3[j]);
					std::thread th4(&MD::Simulation, run4, rho[d], T[t], n[i], A4[j]);

					th1.join(); th2.join(); th3.join(); th4.join();
					delete run1, run2, run3, run4;

					++num;
				}
			}
		}
	}
	// This could be performed at the end of the T-loop above
	// Check stability before implementing
	std::vector<double> a;
	a.reserve(A1.size() + A2.size() + A3.size() + A4.size());
	a.reserve(A1.size() + A2.size());
	a.insert(a.end(), A1.begin(), A1.end());
	a.insert(a.end(), A2.begin(), A2.end());
	a.insert(a.end(), A3.begin(), A3.end());
	a.insert(a.end(), A4.begin(), A4.end());

	// This bit of code has not been tested in its current 2-loop 1-loop format
	for (size_t d = 0; d < rho.size(); d++) {
		for (size_t t = 0; t < T.size(); t++) {
			// Defining in Heap to avoid obj corruption due to multiple optimisation
			Stat_Analysis* data_averaging = new Stat_Analysis(dir_windows, STEPS, PARTICLES, 2.0, 1.5, a);
			for (size_t i = 0; i < n.size(); i++) {
				data_averaging->StaticDataProcessing(n[i]);
			}
			delete data_averaging;
		}
	}

}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
	/*
	* Produces an equally spaced vector of N increments
	* in the inclusive range of [a, b]
	*/
	double h = (b - a) / static_cast<double>(N - 1);
	std::vector<double> xs(N);
	std::vector<double>::iterator x;
	double val;
	for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
		*x = val;
	}
	return xs;
}