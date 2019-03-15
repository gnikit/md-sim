#include <iostream>
#include "phase_transition.h"

#define STEPS_PER_COMPRESSION 5000
#define PARTICLES 500   // 5**3*4
typedef std::vector<double> vec1d;

/* Linux working directory */
std::string dir_linux = ".";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
  vec1d temperature_array = LinearSpacedArray(0.0045, 0.01, 1);
  /* 
   * Adjust the DENSITY, FINAL_DENSITY for the fluid->solid transition and then
   * re-adjust it for the solid->fluid. Otherwise the simulations will sample
   * a lot of unneeded densities.
   */
  for (const auto& i : temperature_array) {
    phase_transition* run1 = new phase_transition(dir_linux, STEPS_PER_COMPRESSION, true, 500, 4, "FCC", false, 1000);
    run1->crystallisation(0.1, 0.2, 0.001, 0.005, 12, 0, "GCM");
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