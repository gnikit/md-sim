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
  size_t steps = 1000;
  MD run(steps, {5, 5, 5}, "SC");
  run.set_visualisation_flag(true);
  run.simulation("3D_view_", 0.5, 0.5, 8, 0.5, "BIP");
}