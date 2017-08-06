#include "MD.h"
//#include "stat_analysis.h"
//#include <memory>
//#include <thread>



int main() {
  size_t steps = 3000;
  size_t num = 1;
  char full = 'y';
  std::string dir = "../../Archives of Data/tests/4^3/";
  std::vector<int> p = { 6/*, 8*/};
  std::vector<double> A = { 0/*, 0.75, 1, 1.25, 1.5 */};
  double density = 0.5;

  MD run(dir, density, steps);
  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < A.size(); j++) {
      srand(time(NULL));
      std::cout << "p: " << p[i] << " A: " << A[j] << " run num: " << num << std::endl;
      run.Simulation(p[i], A[j]);
      ++num;
    }
  }
  //system("pause");
}