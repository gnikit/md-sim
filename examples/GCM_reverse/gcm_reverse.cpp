#include <iostream>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "phase_transition.h"

int main() {
  std::string dir = "./";  // Places files in the execution dir
  size_t compress_steps = 5000;

  /*
   * This method effectively replaces timestep iterations with
   * density iterations by equating the compression timestep to 1
   *
   */
  phase_transition run(dir, compress_steps, true, 500, 7, "FCC", false, 500);

  run.run_backwards("compress_", 0.05, 0.15, 0.05, 0.003, 0, 0, "GCM");
}
