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

  options_type opts;
  opts.simulation_type = "gcm_bcc_";
  opts.potential_type = "GaussianCoreModel";
  opts.lattice = "BCC";
  opts.iterative_method = "VerletAlgorithm";
  opts.particles = {8, 8, 8};
  opts.steps = 100;
  opts.density = 0.5;  //! loop
  opts.target_temperature = 0.0001;
  opts.rdf_options.rdf_bins = 500;
  opts.rdf_options.rdf_wait = 0;

  std::vector<double> rho = helper_functions::linspace(0.05, 0.20, 4);
  for (auto const &i : rho) {
    opts.density = i;
    MD run(opts);
    run.simulation();
  }
}