#include "MD.h"
#include <thread>
#include <string>

//void MultiThread(MD class_object, std::vector<long double> n, std::vector<int> A){
//  for (size_t i = 0; i < n.size(); i++){
//    for (size_t j = 0; j < A.size(); j++){
//      std::cout << "Thread2: p: " << n[i] << " A: " << A[j] << std::endl;
//      class_object.Simulation(n[i], A[j]);
//    }
//  }
//}
//#include "stat_analysis.h"

int main() {
  size_t steps = 10000;
  size_t num = 1;
  double density = 0.5;
  //std::string dir = "/home/gn/Code/C++/MD simulation/Archives of Data/tests/10^3/";
  std::string dir = "../../Archives of Data/tests/10^3/";
  //std::string file = "../../Archives of Data/tests/10^3/Density 0.5/Isothermal~step 10000/";
  std::vector<int> p = { 6 , 8};
  std::vector<int> p2 = { 10, 12 };
  std::vector<double> A = { 0.0/*, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 1.00, 1.05,
                            1.1, 1.2, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00*/ };
  
  system("pause");
  //MD run(dir, density, steps);
  //MD *run2 = new MD(dir, density, steps);
  //
  //for (size_t i = 0; i < p.size(); i++) {
  //  for (size_t j = 0; j < A.size(); j++) {
  //    std::cout << "p: " << p[i] << " A: " << A[j] << " run num: " << num << std::endl;
  //    std::thread th(&MD::Simulation, run2, p2[i], A[j]);
  //    run.Simulation(p[i], A[j]);
  //    th.join();
  //    ++num;
  //  }
  //}
  //delete run2;
  // system("pause");
}
