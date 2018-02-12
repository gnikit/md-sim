#include "MD.h"
#include "stat_analysis.h"
#include <thread>
#include <string>

#define DENSITY 0.5
#define STEPS 10000
#define DIR "../../Archives of Data/tests/10^3/"
#define FILE "../../Archives of Data/tests/10^3/Density 0.5/Isothermal~step 10000/"

//void MultiThread(MD class_object, std::vector<long double> n, std::vector<int> A){
//  for (size_t i = 0; i < n.size(); i++){
//    for (size_t j = 0; j < A.size(); j++){
//      std::cout << "Thread2: p: " << n[i] << " A: " << A[j] << std::endl;
//      class_object.Simulation(n[i], A[j]);
//    }
//  }
//}

int main() {
  size_t steps = 10000;
  size_t num = 1;
  double density = 0.5;
  //std::string dir = "/home/gn/Code/C++/MD simulation/Archives of Data/tests/10^3/";
  std::vector<int> p = { 6, 8, 10, 12 };
  std::vector<double> A1 = { 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85 };
  std::vector<double> A2 = { 0.90, 1.00, 1.05, 1.1, 1.2, 1.25, 1.50 };
  std::vector<double> A3 = { 1.75, 2.00, 2.25, 2.50, 2.75, 3.0, 3.5 };

  MD run(DIR, density, steps);
  MD *run2 = new MD(DIR, density, steps);
  MD *run3 = new MD(DIR, density, steps);

  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < A1.size(); j++) {
      std::cout << "p: " << p[i] << " A: " << A1[j] << " run num: " << num << std::endl;
      std::thread th(&MD::Simulation, run2, p[i], A2[j]);
      std::thread th2(&MD::Simulation, run3, p[i], A3[j]);
      run.Simulation(p[i], A1[j]);
      th.join();
      th2.join();
      ++num;
    }
  }
  delete run2, run3;
  std::vector<double> a;
  a.reserve(A1.size() + A2.size() + A3.size());
  a.insert(a.end(), A1.begin(), A1.end());
  a.insert(a.end(), A2.begin(), A2.end());
  a.insert(a.end(), A3.begin(), A3.end());
  a.insert(a.begin(), 0);
  a.push_back(4.);

  Stat_Analysis test(FILE, a);
  for (size_t i = 0; i < p.size(); i++) {
    test.StaticDataProcessing(p[i]);
  }
  //system("pause");
}