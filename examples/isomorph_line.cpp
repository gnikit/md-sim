#include <iostream>
#include <string>
#include <thread>
#include "MD.h"
#include "isomorph.h"

#define STEPS 15000
#define PARTICLES 1000
typedef std::vector<double> vec1d;

std::string dir_linux = "/home/gn/Desktop/test_data/isomorph";

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

int main() {
  /* Potential power strength */
  size_t n = 8;
  /* Generate Temperature vector for isomorph */
  vec1d T_iso = LinearSpacedArray(0.5, 2, 5);
  /* Empty containers for density and a_par */
  vec1d rho_iso, A_iso;
  vec1d rho_iso_h, A_iso_h;
  vec1d rho_iso_k, A_iso_k;
  vec1d rho_iso_l, A_iso_l;
  /* Generate Density and A isomorph vectors */
  Isomorph isomorph_line(0.5, 0.5, 0.5, T_iso);
  std::tie(rho_iso, A_iso) = isomorph_line.GenLine(n);
  Isomorph isomorph_line_h(0.5, 0.5, 0.75, T_iso);
  std::tie(rho_iso_h, A_iso_h) = isomorph_line_h.GenLine(n);
  Isomorph isomorph_line_k(0.5, 0.5, 1.25, T_iso);
  std::tie(rho_iso_k, A_iso_k) = isomorph_line_k.GenLine(n);
  Isomorph isomorph_linr_l(0.5, 0.5, 2.00, T_iso);
  std::tie(rho_iso_l, A_iso_l) = isomorph_linr_l.GenLine(n);
  /*
   * This is an isomorph line run
   * Simulates the fluid along the line
   */
  for (size_t i = 0; i < T_iso.size(); ++i) {
    MD run(dir_linux, STEPS, false, 500, 10, false, 2000);
    run.Simulation(rho_iso[i], T_iso[i], n, A_iso[i]);
  }
}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
  /*
   * Produces an equally spaced vector of N increments
   * in the inclusive range of [a, b]
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