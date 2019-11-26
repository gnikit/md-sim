#include <string>
#include "MD.h"

/*
 * A demonstration of how the MD class should be used.
 * It creates 3 files in the directory.
 * One containing all the state properties of the fluid,
 * along with the Mean Square Displacement and the
 * Velocity Autocorrelation Function.
 * A second file containing the Radial Distribution Function
 * and a third file containing the positions, velocities
 * and accelerations of the particles at the last time step.
 */

int main() {
  std::string dir = ".";
  size_t steps = 5000;
  MD run(dir, steps, false, 500, 10, "SC", false, 2000);
  run.simulation("lj_simple_run_", 0.5, 0.5, NAN, NAN, "LJ");
}