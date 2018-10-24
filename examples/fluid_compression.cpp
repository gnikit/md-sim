#include <string>
#include <thread>
#include "MD.h"

#define STEPS 1000
#define PARTICLES 1000
typedef std::vector<double> vec1d;

/* Linux working directory */
std::string dir_linux = "/home/gn/Code/MD-simulation/examples/example_data";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
  vec1d temperature_array = LinearSpacedArray(0.001, 0.01, 10);
  for (const auto& i : temperature_array) {
    MD* run1 = new MD(dir_linux, STEPS, true);  // fluid compression set to true
    MD* run2 = new MD(dir_linux, STEPS, true);
    MD* run3 = new MD(dir_linux, STEPS, true);
    MD* run4 = new MD(dir_linux, STEPS, true);

    std::thread th1(&MD::Simulation, run1, 0.05, i, 12, 0.25);
    std::thread th2(&MD::Simulation, run2, 0.05, i, 12, 0.5);
    std::thread th3(&MD::Simulation, run3, 0.05, i, 12, 0.1);
    std::thread th4(&MD::Simulation, run4, 0.05, i, 12, 1.25);

    th1.join();
    th2.join();
    th3.join();
    th4.join();

    delete run1;
    delete run2;
    delete run3;
    delete run4;
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