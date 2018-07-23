#include <string>
#include <thread>
#include "MD.h"
#include "isomorph.h"
#include "stat_analysis.h"

#define STEPS 1000      //10000
#define PARTICLES 1000  //1000
typedef std::vector<double> vec1d;

/* Linux working directory */
std::string dir_linux = "/home/gn/Desktop/test_data/delete/";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
  vec1d temp = LinearSpacedArray(0.001, 0.01, 10);
  for (size_t i = 0; i < temp.size(); i++) {
  	MD* run1 = new MD(dir_linux, STEPS, false);  // fluid compression set to false
    MD* run2 = new MD(dir_linux, STEPS, false);
    MD* run3 = new MD(dir_linux, STEPS, false);
    MD* run4 = new MD(dir_linux, STEPS, false);

    std::thread th1(&MD::Simulation, run1, 0.05, 0.5, 12, 0.0);
    // std::thread th2(&MD::Simulation, run2, 0.05, 0.5, 12, 0.5);
    // std::thread th3(&MD::Simulation, run3, 0.05, 0.5, 12, 0.1);
    // std::thread th4(&MD::Simulation, run4, 0.05, 0.5, 12, 1.5);

    th1.join(); /* th2.join(); th3.join(); th4.join(); */
    delete run1, run2, run3, run4;
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