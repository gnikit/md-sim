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
  phase_transition run(compress_steps, 7, "FCC");
  run.set_compression_flag(true);

  std::vector<double> temperatures = {
      4.602932914279152937e-05, 6.461042397802770065e-05,
      6.454528429411781482e-05, 7.688925439500345071e-05,
      8.923322449588821925e-05, 1.078306042521020788e-04};

  for (size_t t = 0; t < temperatures.size(); ++t) {
    run.crystallisation("compress_fwd_" + std::to_string(t), 3.40e-02, 4.2e-02,
                        double(abs(3.40e-02 - 4.2e-02) / 8),
                        temperatures[t], NAN, NAN, "GCM");
  }
}
