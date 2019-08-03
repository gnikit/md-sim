#include <iostream>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "phase_transition.h"

int main() {
  std::string dir = "/home/gn/Desktop/md_data/compress_gcm";
  size_t compress_steps = 100;
  //todo: create a range density vector
  //todo: loop through it one timestep at a time
  /*
    * This method effectively replaces timestep iterations with 
    * density iterations by equating the compression timestep to 1
    * 
   */
  phase_transition run(dir, compress_steps, true, 10, 6, "FCC", false, 0);
  std::vector<double> temperatures = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006};
  for (int t = 0; t < temperatures.size(); ++t) {
    run.crystallisation(0.05, 0.8, 0.001, temperatures[t], 0, 0, "GCM");
  }
}
