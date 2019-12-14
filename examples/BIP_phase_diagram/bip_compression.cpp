#include <iostream>
#include "helper_functions.h"
#include "phase_transition.h"

#define STEPS_PER_COMPRESSION 5000
typedef std::vector<double> vec1d;

int main() {
  vec1d temperature_array = helper_functions::linspace(0.006, 0.01, 1);
  /*
   * Adjust the DENSITY, FINAL_DENSITY for the fluid->solid transition and then
   * re-adjust it for the solid->fluid. Otherwise the simulations will sample
   * a lot of unneeded densities.
   */

  for (const auto& i : temperature_array) {
    phase_transition run(STEPS_PER_COMPRESSION, {10, 10, 10}, "SC");
    run.set_compression_flag(true);
    run.crystallisation("compress_fwd_", 0.05, 0.1, 0.01, i, 12, 0, "BIP");
  }
}
