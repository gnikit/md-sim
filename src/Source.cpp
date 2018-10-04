#include <string>
#include <thread>
#include "MD.h"

#define STEPS 10000     //10000

/* Linux working directory */
std::string dir = "/home/gn/Desktop/test_data/sample/";

int main(int argc, char const *argv[]) {
  MD run(dir, STEPS);
  run.Simulation(0.5, 0.5, 6, 0.5);
}
