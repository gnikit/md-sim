#include "MD.h"
#include "isomorph.h"
#include "stat_analysis.h"
#include <thread>
#include <string>

#define STEPS 10000
#define PARTICLES 1000
typedef std::vector<double> vec1d;

/* Windows working directory for Archives of Data */
std::string dir_windows = "C:/Code/C++/MD simulation/Archives of Data/testing/";
/* Linux working directory */
std::string dir_linux = "/home/gn/Desktop/test_data/";
std::string dir_crystal = "/home/gn/Desktop/crystallisation_data/";
std::string dir_crystal_win = "C:/Code/C++/MD simulation/Archives of Data/crystallisation/gaussian/";
/* Working directory of the cx1 cluster */
std::string dir = "";

void MakeDataBase();

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {

	/* Potential power strength */
	size_t n = 8;
	/* Generate Temperature vector for isomorph */
	vec1d T_iso = LinearSpacedArray(0.5, 3, 5);
	/* Empty containers for density and a_par */
	vec1d rho_iso, A_iso;
	vec1d rho_iso_h, A_iso_h;
	vec1d rho_iso_k, A_iso_k;
	vec1d rho_iso_l, A_iso_l;
	/* Generate Density and A isomorph vectors */
	Isomorph isomorph_line(0.5, 0.5, 1, T_iso);
	std::tie(rho_iso, A_iso) = isomorph_line.GenLine(n);

	Isomorph isomorph_line_h(0.5, 0.5, 0.75, T_iso);
	std::tie(rho_iso_h, A_iso_h) = isomorph_line_h.GenLine(n);

	Isomorph isomorph_line_k(0.5, 0.5, 1.25, T_iso);
	std::tie(rho_iso_k, A_iso_k) = isomorph_line_k.GenLine(n);

	Isomorph isomorph_linr_l(0.5, 0.5, 2.00, T_iso);
	std::tie(rho_iso_l, A_iso_l) = isomorph_linr_l.GenLine(n);

	/* Simulation Examples */
	//MD run(dir_linux, STEPS);
	//run.Simulation(0.5, 0.5, 8, 0.5);
	/*
	* This is an isomorph line run
	* Simulates the fluid along the line
	*/
	/* for (size_t i = 0; i < T_iso.size(); i++) {
	std::cout << "T: " << T_iso[i] << " rho: " << rho_iso[i] << " A: " << A_iso[i] << std::endl;
	MD run(dir_windows, STEPS);
	run.Simulation(rho_iso[i], T_iso[i], n, A_iso[i]);
	}*/
	/* Individual Runs */
	MD* run1 = new MD(dir_windows, STEPS);
	MD* run2 = new MD(dir_windows, STEPS);
	MD* run3 = new MD(dir_windows, STEPS);
	MD* run4 = new MD(dir_windows, STEPS);
	vec1d a = { 0, 0.25 };

	for (size_t i = 0; i < a.size(); i++) {
		std::thread th1(&MD::Simulation, run1, 0.5, 1, 6, a[i]);
		th1.join();
	}
	//std::thread th2(&MD::Simulation, run2, 0.5, 0.5, 12,  0.75);
	// std::thread th3(&MD::Simulation, run3, 0.6, 10, 10, 0.5);
	// std::thread th4(&MD::Simulation, run4, 0.6, 10, 12, 0.5);

	// th2.join();// th3.join(); th4.join();
	delete run1, run2, run3, run4;

	/*-----------------------------------------------*/

	//system("pause");
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