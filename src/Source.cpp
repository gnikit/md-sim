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
  std::string dir = "../../Archives of Data/";  // Current Working Directory
  std::vector<int> p = { 6, 8, 10, 12 };
  double a_value = 0.5 * std::pow(0.5, (2 / 12));
  // std::vector<double> A1 = { 0, 0.25, 0.50, 0.75, a_value, 1.00, 1.25, 1.50, 1.75, 2.00, 3, 4.00 };
  // std::vector<double> A1 = { 0.0/* , 0.25, 0.50, 0.65, 0.70, 0.75 */};
  // std::vector<double> A2 = { 0.5/* , 0.75, 0.90, 0.95, 1.00, 1.05 */};
  // std::vector<double> A3 = { 1.0/* , 1.25, 1.25, 1.50, 1.75, 2.00 */};
  // std::vector<double> A4 = { 2.0/* , 2.50, 2.75, 3.00, 3.50, 4.00 */};
  std::vector<double> A1 = LinearSpacedArray(0,1,5);
  std::vector<double> A2 = LinearSpacedArray(1,2,5);
  std::vector<double> A3 = LinearSpacedArray(2,4,5);
  std::vector<double> A4 = LinearSpacedArray(5,10,5);
  std::vector<double> A5 = LinearSpacedArray(0,1,5);
  std::vector<double> A6 = LinearSpacedArray(1,2,5);
  std::vector<double> A7 = LinearSpacedArray(2,4,5);
  std::vector<double> A8 = LinearSpacedArray(5,10,5);

  MD run(dir, STEPS);
  MD* run2 = new MD(dir, STEPS);
  MD* run3 = new MD(dir, STEPS);
  MD* run4 = new MD(dir, STEPS);
  MD* run5 = new MD(dir, STEPS);
  MD* run6 = new MD(dir, STEPS);
  MD* run7 = new MD(dir, STEPS);
  MD* run8 = new MD(dir, STEPS);

  for (size_t i = 0; i < p.size(); i++) {
    for (size_t j = 0; j < A1.size(); j++) {
      std::cout << "p: " << p[i] << " A: " << A1[j] << " run num: " << num << std::endl;
      // std::thread th2(&MD::Simulation, run2, p[i], A2[j]);
      // std::thread th3(&MD::Simulation, run3, p[i], A3[j]);
      // std::thread th4(&MD::Simulation, run4, p[i], A4[j]);
      // std::thread th5(&MD::Simulation, run5, p[i], A5[j]);
      // std::thread th6(&MD::Simulation, run6, p[i], A6[j]);
      // std::thread th7(&MD::Simulation, run7, p[i], A7[j]);
      // std::thread th8(&MD::Simulation, run8, p[i], A8[j]);
      run.Simulation(0.5, 0.5, 8, 0.5);

      // th2.join();
      // th3.join();
      // th4.join();
      // th5.join();
      // th6.join();
      // th7.join();
      // th8.join();
      ++num;
    }
  }
  delete run2, run3, run4, run5, run6, run7, run8;
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
  //system("pause");
}