#include <string>
#include <thread>
#include "MD.h"
#include "isomorph.h"
#include "stat_analysis.h"

#define STEPS 10000     //10000
#define PARTICLES 1000  //1000
typedef std::vector<double> vec1d;

/* Linux working directory */
std::string dir_linux = "/home/gn/Desktop/test_data/sample/";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
  // size_t num = 1;
  // std::vector<size_t> n = {6, 8, 10, 12};
  // std::vector<double> rho = {/* 0.2, */ 0.5, 1.0};
  // std::vector<double> T = {0.5, 1.0, 2.0};
  // std::vector<double> A1 = LinearSpacedArray(0, 1, 5);
  // std::vector<double> A2 = LinearSpacedArray(1.25, 2.25, 5);
  // // std::vector<double> A3 = LinearSpacedArray(2.50, 4.50, 5);
  // std::vector<double> A4 = {0.70, 0.80, 0.90, 0.95, 1.10};

  // for (size_t d = 0; d < rho.size(); d++) {
  //   for (size_t t = 0; t < T.size(); t++) {
  //     for (size_t i = 0; i < n.size(); i++) {
  //       for (size_t j = 0; j < A1.size(); j++) {
  //         std::cout << "rho: " << rho[d] << " T: " << T[t] << "n: " << n[i]
  //                   << " A: " << A1[j] << " run num: " << num << std::endl;

  //         MD* run1 = new MD(dir_linux, STEPS);
  //         MD* run2 = new MD(dir_linux, STEPS);
  //         MD* run3 = new MD(dir_linux, STEPS);
  //         MD* run4 = new MD(dir_linux, STEPS);

  //         std::thread th1(&MD::Simulation, run1, rho[d], T[t], n[i], A1[j]);
  //         std::thread th2(&MD::Simulation, run2, rho[d], T[t], n[i], A2[j]);
  //         // std::thread th3(&MD::Simulation, run3, rho[d], T[t], n[i], A3[j]);
  //         std::thread th4(&MD::Simulation, run4, rho[d], T[t], n[i], A4[j]);

  //         th1.join();
  //         th2.join();
  //         // th3.join();
  //         th4.join();
  //         delete run1, run2, run3, run4;

  //         ++num;
  //       }
  //     }
  //   }
  // }
  MD run1(dir_linux, STEPS);
  run1.Simulation(0.5, 0.5, 6, 0.5);

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
