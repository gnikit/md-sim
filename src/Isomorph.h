#pragma once
#include <tuple>
#include <vector>


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
  Isomorph(double RHO, double T, double Ar, vec1d T_in);

  double getRho(double rho1, double T1, double T2, size_t n);

  double getA(double a1, double rho1, double rho2, size_t n);
  std::tuple<vec1d, vec1d> GenLine(size_t n);
};