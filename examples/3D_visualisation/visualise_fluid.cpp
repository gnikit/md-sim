#include <string>
#include "MD.h"

/*
 * This example creates the normal data files that can be used
 * to statistically analyse the fluid. In addition, to those files,
 * 3 files containing the positions of the particles over time are also
 * created. The script visualise_fluid.py in the tools directory can be
 * used to produce a 3D scatter plot animation of the fluid.
 */

int main() {
  std::string dir = ".";
  size_t steps = 5000;
  MD run(dir, steps, false, 100, 10, true, 0);
  run.Simulation(0.5, 0.5, 8, 0.5, "BIP");
}