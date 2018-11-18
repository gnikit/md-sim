#include <thread>
#include "MD.h"
// BUG: potential bug, including 2 different math libraries if intel compiler
#include <math.h>
#include <iostream>
#define PARTICLES 1000
typedef std::vector<double> vec1d;

/* Linux working directory */
std::string dir_linux = "/home/gn/Desktop/test_data/gaussian";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

// TODO: Insert the individual number of steps per compression and the initial
// and final compression limits, along with increments. Then let the simulation
// run
int main() {
  vec1d temperature_array = LinearSpacedArray(0.001, 0.01, 10);

  unsigned int COMPRESS_EVERY = 1000;
  double rho = 0.05;
  double final_rho = 0.8;
  double rho_increment = 0.001;
  unsigned int compression_num = ceil((final_rho - rho) / rho_increment);
  unsigned int STEPS = compression_num * COMPRESS_EVERY;

  std::cout << "Initial density: " << rho << "\nFinal density: " << final_rho
            << "\nDensity increment: " << rho_increment
            << "\nTotal steps: " << STEPS
            << "\nCompress every: " << COMPRESS_EVERY << std::endl;

  for (const auto& i : temperature_array) {
    MD* run1 = new MD(dir_linux, STEPS, true);  // fluid compression set to true
    run1->steps_per_compress = COMPRESS_EVERY;
    run1->Simulation(rho, i, 12, 0);
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