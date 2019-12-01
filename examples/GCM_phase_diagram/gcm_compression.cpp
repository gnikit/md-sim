#include <iostream>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "phase_transition.h"

int main() {
  size_t compress_steps = 5000;

  /*
   * This method effectively replaces timestep iterations with
   * density iterations by equating the compression timestep to 1
   *
   */
  phase_transition run(compress_steps, 5, "FCC");
  run.set_compression_flag(true);

  std::vector<double> temperatures = {0.0033, 0.0035};
  for (size_t t = 0; t < temperatures.size(); ++t) {
    run.crystallisation("compress_fwd_", 0.05, 0.20, 0.025, temperatures[t],
                        NAN, NAN, "GCM");
  }
}
