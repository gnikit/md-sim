#include <iostream>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "phase_transition.h"

int main() {
  size_t compress_steps = 1000;

  /*
   * This method effectively replaces timestep iterations with
   * density iterations by equating the compression timestep to 1
   *
   */
  phase_transition run(compress_steps, 7, "FCC");

  double rho_start = 1.0e-01;
  double rho_end = 2.0e-01;
  std::vector<double> t2 = helper_functions::linspace(5.0e-04, 3.0e-03, 20);

  for (size_t t = 0; t < t2.size(); ++t) {
    run.crystallisation("compress_fwd_" + std::to_string(t), rho_start, rho_end,
                        (rho_end - rho_start) / 20, t2[t], NAN, NAN, "GCM");
  }
}
