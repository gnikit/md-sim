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
  options_type opts;
  opts.steps = 1000;
  opts.particles = {5, 5, 5};
  opts.lattice = "SC";
  opts.potential_type = "BoundedInversePower";
  opts.density = 0.5;
  opts.target_temperature = 0.5;
  opts.power = 8;
  opts.io_options.simulation_name = "3D_view_";
  
  MD run(opts);
  run.simulation();
}