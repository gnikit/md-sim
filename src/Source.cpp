#include "MD.h"
#include "stat_analysis.h"
#include <thread>
#include <string>
#include <tuple>

#define STEPS 10000
typedef std::vector<double> vec1d;

// TODO: Move to another file?
class Isomorph {
  typedef std::vector<double> vec1d;
  vec1d _T, _RHO, _A; // T points
  double _rho_r;      // Reference density
  double _T_r;        // Reference temperature
  double _A_r;        // Reference A
  /*----------------------------------*/
  double _rho_out;    // Output density
  double _T_out;      // Output temperature
  double _A_out;      // Output A
  public:
  Isomorph(double RHO, double T, double Ar, vec1d T_in) {
    /*
    * Takes as arguments a, Reference density, temperature and A parameter
    * T_in is the input Temperature vector, where the isomorph will be placed
    */
    _rho_r = RHO;
    _T_r = T;
    _A_r = Ar;
    _T = T_in;
  }

  double getRho(double rho1, double T1, double T2, size_t n) {
    double rho2 = rho1 * pow((T2 / T1), (3.0 / n));
    return rho2;
  }

  double getA(double a1, double rho1, double rho2, size_t n) {
    double a2 = a1 * pow((rho1 / rho2), (1.0 / 3.0));
    return a2;
  }
  std::tuple<vec1d, vec1d> GenLine(size_t n) {
    /*
    * This method returns a tuple of Isomorphic points
    * for the density and A parameter in the form of a vector.
    */
    for (size_t i = 0; i < _T.size(); i++) {
      _T_out = _T[i];   // reduntant step
      _rho_out = getRho(_rho_r, _T_r, _T_out, n);
      _A_out = getA(_A_r, _rho_r, _rho_out, n);
      _RHO.push_back(_rho_out);
      _A.push_back(_A_out);
    }
    return std::make_tuple(_RHO, _A);
  }
};


std::vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
  /*
  * Produces an equally spaced vector of N increments
  * in the inclusive range of (a, b)
  */
  double h = (b - a) / static_cast<double>(N - 1);
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
  /* Windows working directory for Archives of Data */
  std::string dir_windows = "C:/Code/C++/MD simulation/Archives of Data/"; 
  /* Working directory of the cx1 cluster */
  std::string dir = "/home/gn/test_data/";
  /* Potential power strength */
  size_t n = 8;
  /* Generate Temperature vector for isomorph */
  vec1d T_iso = LinearSpacedArray(0.5, 3, 5);
  /* Empty containers for density and a_par */
  vec1d rho_iso, A_iso;
  vec1d rho_iso_h, A_iso_h;
  vec1d rho_iso_k, A_iso_k;
  vec1d rho_iso_l, A_iso_l;
  /* Generate Density and A isomorph vectors */
  Isomorph isomorph_line(0.5, 0.5, 1, T_iso);
  std::tie(rho_iso, A_iso) = isomorph_line.GenLine(n);

  Isomorph isomorph_line_h(0.5, 0.5, 0.75, T_iso);
  std::tie(rho_iso_h, A_iso_h) = isomorph_line_h.GenLine(n);

  Isomorph isomorph_line_k(0.5, 0.5, 1.25, T_iso);
  std::tie(rho_iso_k, A_iso_k) = isomorph_line_k.GenLine(n);

  Isomorph isomorph_linr_l(0.5, 0.5, 2.00, T_iso);
  std::tie(rho_iso_l, A_iso_l) = isomorph_linr_l.GenLine(n);
  MD run(dir, STEPS);
  run.Simulation(0.5, 0.5, 8, 0.5);
  /*
  * This is an isomorph line run
  * Simulates the fluid along the line
  */
 /* for (size_t i = 0; i < T_iso.size(); i++) {
    std::cout << "T: " << T_iso[i] << " rho: " << rho_iso[i] << " A: " << A_iso[i] << std::endl;
    MD run(dir_windows, STEPS);
    run.Simulation(rho_iso[i], T_iso[i], n, A_iso[i]);
  }*/
  /* Individual Run */
  //MD* run1 = new MD(dir_windows, STEPS);
  //MD* run2 = new MD(dir_windows, STEPS);
  //MD* run3 = new MD(dir_windows, STEPS);
  //MD* run4 = new MD(dir_windows, STEPS);

  //std::thread th1(&MD::Simulation, run1, 1, 1, 6,  0.5);
  //std::thread th2(&MD::Simulation, run2, 1, 1, 8,  0.5);
  //std::thread th3(&MD::Simulation, run3, 1, 1, 10, 0.5);
  //std::thread th4(&MD::Simulation, run4, 1, 1, 12, 0.5);

  //th1.join(); th2.join(); th3.join(); th4.join();
  //delete run1, run2, run3, run4;

  /*-----------------------------------------------*/
  //std::vector<size_t> n = { 6, 8, 10, 12 };
  //std::vector<double> rho = { 0.5, 1.0, 1.5, 2.0 }; 
  //std::vector<double> T = { 0.5, 1.0, 1.5, 2.0  };  
  //std::vector<double> A1 = LinearSpacedArray(0, 1, 5);
  //std::vector<double> A2 = LinearSpacedArray(1.25, 2.25, 5);
  //std::vector<double> A3 = LinearSpacedArray(2.50, 4.50, 5);  
  //std::vector<double> A4 = LinearSpacedArray(0.6, 1.2, 5);       
  //
  //for (size_t d = 0; d < rho.size(); d++) {
  //  for (size_t t = 0; t < T.size(); t++) {
  //    for (size_t i = 0; i < n.size(); i++) {
  //      for (size_t j = 0; j < A1.size(); j++) {
  //        std::cout << "p: " << n[i] << " A: " << A1[j] << " run num: " << num << std::endl;
  //
  //        MD* run1 = new MD(dir, STEPS);
  //        MD* run2 = new MD(dir, STEPS);
  //        MD* run3 = new MD(dir, STEPS);
  //        MD* run4 = new MD(dir, STEPS);
  //
  //        std::thread th1(&MD::Simulation, run1, rho[d], T[t], n[i], A1[j]);
  //        std::thread th2(&MD::Simulation, run2, rho[d], T[t], n[i], A2[j]);
  //        std::thread th3(&MD::Simulation, run3, rho[d], T[t], n[i], A3[j]);
  //        std::thread th4(&MD::Simulation, run4, rho[d], T[t], n[i], A4[j]);
  //
  //        th1.join(); th2.join(); th3.join(); th4.join();
  //        delete run1, run2, run3, run4;
  //
  //        ++num;
  //      }
  //    }
  //  }
  //}
  //std::vector<double> a;
  //a.reserve(A1.size() + A2.size() /*+ A3.size() + A4.size()*/);
  //a.reserve(A1.size() + A2.size());
  //a.insert(a.end(), A1.begin(), A1.end());
  //a.insert(a.end(), A2.begin(), A2.end());
  ////a.insert(a.end(), A3.begin(), A3.end());
  ////a.insert(a.end(), A4.begin(), A4.end());
  //
  //Stat_Analysis test(dir_windows, STEPS, 1000, 2.0, 1.5, a);
  // for (size_t i = 0; i < n.size(); i++) {
  //   test.StaticDataProcessing(n[i]);
  // }
  //system("pause");
}