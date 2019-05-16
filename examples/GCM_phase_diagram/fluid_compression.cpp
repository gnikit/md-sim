#include <iostream>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "phase_transition.h"

int main() {
  std::string dir = "/home/gn/Desktop/md_data/compress_gcm";
  size_t compress_steps = 5000;

  phase_transition run(dir, compress_steps, true, 500, 6, "FCC", false, 500);
  std::vector<double> temperatures = {0.001, 0.003, 0.0033, 0.0035, 0.0038};
  for (int t = 0; t < temperatures.size(); ++t) {
    run.crystallisation(0.05, 0.30, 0.025, temperatures[t], 0, 0, "GCM");
  }
}
