#include "MD.h"
#include "stat_analysis.h"




int main() {
  size_t steps = 10000;
  size_t num = 1;
  std::string dir = "../../Archives of Data/tests/10^3/";
  std::string file = "../../Archives of Data/tests/10^3/Density 0.5/Isothermal~step 10000/";
  typedef std::vector<long double> vec1d;
  std::vector<int> p = { 6/*, 8, 10, 12*/ };
  vec1d A = { 0.75, 1, 1.25, 1.50, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 1.00, 1.05,
               1.1, 1.2, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00 };
  long double density = 0.5;

  MD run(dir, density, steps);
  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < A.size(); j++) {
      srand(time(NULL));
      std::cout << "p: " << p[i] << " A: " << A[j] << " run num: " << num << std::endl;
      run.Simulation(p[i], A[j]);
      ++num;
    }
  }
  //Stat_Analysis test(file, a);
  //for (size_t i = 0; i < p.size(); i++) {
  //  test.StaticDataProcessing(p[i]);
  //}
}
