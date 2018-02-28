#include "MD.h"
// #include "stat_analysis.h"
#include <thread>
#include <string>

#define DENSITY 0.5
#define STEPS 10000

//void MultiThread(MD class_object, std::vector<long double> n, std::vector<int> A){
//  for (size_t i = 0; i < n.size(); i++){
//    for (size_t j = 0; j < A.size(); j++){
//      std::cout << "Thread2: p: " << n[i] << " A: " << A[j] << std::endl;
//      class_object.Simulation(n[i], A[j]);
//    }
//  }
//}
double getRho2(double rho1, double T1, double T2, size_t n){
  double rho2 = rho1 * std::pow((T2/T1), (3.0/n));
  return rho2;
}

double getA2(double a1, double rho1, double rho2, size_t n){
  double a2 = a1 * std::pow((rho1/rho2), (1.0/n));
  return a2;
}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N)
    {
        double h = (b - a) / static_cast<double>(N-1);
        std::vector<double> xs(N);
        std::vector<double>::iterator x;
        double val;
        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
            *x = val;
        }
        return xs;
    }

int main() {
  size_t num = 1;
  std::string dir_windows = "C:/Code/C++/MD simulation/Archives of Data/";  // Current Working Directory
  std::string dir = "";   // Working directory of the cluster
  std::vector<size_t> n = { 6/*, 8, 10, 12*/ };
  std::vector<double> rho = { 0.5/*, 1.0, 1.5, 2.0, 2.5*/ };
  std::vector<double> T = { 0.5/*, 1.0, 1.5, 2.0 */};
  //std::vector<double> A1 = { 0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 4.00 };
  std::vector<double> A1 = LinearSpacedArray(0,1,5);
  std::vector<double> A2 = LinearSpacedArray(1.25,2.25,5);
  std::vector<double> A3 = LinearSpacedArray(2.50,4.50,5);
  std::vector<double> A4 = LinearSpacedArray(5,10,5);


  MD run(dir_windows, STEPS);
  MD* run2 = new MD(dir, STEPS);
  MD* run3 = new MD(dir, STEPS);
  MD* run4 = new MD(dir, STEPS);


  for (size_t d = 0; d < rho.size(); d++) {
    for (size_t t = 0; t < T.size(); t++) {
      for (size_t i = 0; i < n.size(); i++) {
        for (size_t j = 0; j < A1.size(); j++) {
          //std::cout << "p: " << n[i] << " A: " << A1[j] << " run num: " << num << std::endl;
          //std::thread th2(&MD::Simulation, run2, rho.at(d), T.at(t), n.at(i), A2.at(j));
          //std::thread th3(&MD::Simulation, run3, rho.at(d), T.at(t), n.at(i), A3.at(j));
          //std::thread th4(&MD::Simulation, run4, rho.at(d), T.at(t), n.at(i), A4.at(j));
          run.Simulation(rho.at(d), T.at(t), n.at(i), A1.at(j));
          
          //th2.join(); th3.join(); th4.join();
          ++num;
        }
      }
    }
  }
  delete run2, run3, run4;
  //std::vector<double> a;
  //a.reserve(A1.size() + A2.size() + A3.size() + A4.size());
  //a.reserve(A1.size() + A2.size());
  //a.insert(a.end(), A1.begin(), A1.end());
  //a.insert(a.end(), A2.begin(), A2.end());
  //a.insert(a.end(), A3.begin(), A3.end());
  //a.insert(a.end(), A4.begin(), A4.end());
  //
  // Stat_Analysis test(dir, A1, STEPS, 1.0, 1000, DENSITY);
  // for (size_t i = 0; i < p.size(); i++) {
  //   test.StaticDataProcessing(p[i]);
  // }
  system("pause");
}