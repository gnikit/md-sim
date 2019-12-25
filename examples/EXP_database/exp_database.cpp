#include <iostream>
#include <string>
#include <thread>
#include "MD.h"
#include "helper_functions.h"

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

#define STEPS 10000

void MakeDataBase();

int main() { MakeDataBase(); }

void MakeDataBase() {
  /*
   * Builds a relatively complete database mapping the rho,T,n,a space
   * of te MD fluid. Performs Also statistical averaging at the end with
   * quantities such as Pc, U, K, E, Pk.
   * NOTE: It takes a long time to complete and generates multiple files!
   */

  size_t num = 1;
  std::vector<double> n = {1.5, 2.0, 2.5, 3.0};
  std::vector<double> rho = {0.20, 0.3, 0.5, 0.8, 1.0};

  std::vector<double> T = {0.5, 1.0, 1.5, 2.0};
  std::vector<double> A1 = helper_functions::linspace(0.0, 1.0, 5);
  std::vector<double> A2 = helper_functions::linspace(1.25, 2.25, 5);
  std::vector<double> A3 = helper_functions::linspace(2.50, 4.50, 5);
  std::vector<double> A4 = {0.20, 0.40, 0.60, 0.80, 0.90};

  for (size_t d = 0; d < rho.size(); ++d) {
    for (size_t t = 0; t < T.size(); ++t) {
      for (size_t i = 0; i < n.size(); ++i) {
        std::cout << " run num: " << num << std::endl;

        MD* run1 = new MD(STEPS, {8, 8, 8}, "SC");
        MD* run2 = new MD(STEPS, {8, 8, 8}, "SC");

        // ? set temperature manually in default machine
        std::thread th1(
            static_cast<void (MD::*)(std::string, double, double, double,
                                     double, std::string)>(&MD::simulation),
            run1, "", rho[d], T[t], n[i], exp(0.25), "Exponential");
        std::thread th2(
            static_cast<void (MD::*)(std::string, double, double, double,
                                     double, std::string)>(&MD::simulation),
            run2, "", rho[d], T[t], n[i], exp(0.75), "Exponential");

        th1.join();
        th2.join();

        delete run1;
        delete run2;

        ++num;
      }
    }
  }
}
