#include <iostream>
#include <string>
#include <thread>
#include "MD.h"

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

#define STEPS 65000  // 10000
#define RDF_EQ 15000

std::string dir_linux = "/home/gn/Desktop/md_data/exp";

void MakeDataBase();

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() { MakeDataBase(); }

void MakeDataBase() {
  /*
   * Builds a relatively complete database mapping the rho,T,n,a space
   * of te MD fluid. Performs Also statistical averaging at the end with
   * quantities such as Pc, U, K, E, Pk.
   * NOTE: It takes a long time to complete and generates multiple files!
   */

  size_t num = 1;
  std::vector<double> n = {1.5, 2.0, 2.5, 3.0 /* 8, 10, 12 */};
  std::vector<double> rho = {0.20, 0.3, 0.5, 0.8, 1.0};
  /* Test temperature smaller than 1 line 0.75 */
  std::vector<double> T = {1 /* 0.5, 1.0, 1.5, 2.0 */};
  std::vector<double> A1 = LinearSpacedArray(0, 1, 5);
  std::vector<double> A2 = LinearSpacedArray(1.25, 2.25, 5);
  std::vector<double> A3 = LinearSpacedArray(2.50, 4.50, 5);
  std::vector<double> A4 = {0.20, 0.40, 0.60, 0.80, 0.90};

  for (size_t d = 0; d < rho.size(); ++d) {
    for (size_t t = 0; t < T.size(); ++t) {
      for (size_t i = 0; i < n.size(); ++i) {
        // for (size_t j = 0; j < A1.size(); ++j) {
        std::cout << " run num: " << num << std::endl;

        MD* run1 = new MD(dir_linux, STEPS, false, 250, 8, "SC", false, RDF_EQ);
        MD* run2 = new MD(dir_linux, STEPS, false, 250, 8, "SC", false, RDF_EQ);
        // MD* run3 = new MD(dir_linux, STEPS, false, 500, 10, false, RDF_EQ);
        // MD* run4 = new MD(dir_linux, STEPS, false, 500, 10, false, RDF_EQ);

        // ? set temperature manually in default machine
        std::thread th1(&MD::Simulation, run1, rho[d], T[t], n[i], exp(0.25),
                        "EXP");
        std::thread th2(&MD::Simulation, run2, rho[d], T[t], n[i], exp(0.75),
                        "EXP");
        // std::thread th3(&MD::Simulation, run3, rho[d], T[t], n[i], A3[j]);
        // std::thread th4(&MD::Simulation, run4, rho[d], T[t], n[i], A4[j]);

        th1.join();
        th2.join();
        // th3.join();
        // th4.join();

        delete run1;
        delete run2;
        // delete run3;
        // delete run4;

        ++num;
        // }
      }
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