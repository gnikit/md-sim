#include <iostream>
#include <string>
#include <thread>
#include "MD.h"
#include "helper_functions.h"
#include "isomorph.h"

#define STEPS 5000
#define PARTICLES 1000
typedef std::vector<double> vec1d;

int main() {
  /* Potential power strength */
  size_t n = 8;
  /* Generate Temperature vector for isomorph */
  vec1d T_iso = helper_functions::linspace(0.5, 2.0, 5);
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
    MD run(STEPS, {8, 8, 8}, "SC");
    run.simulation("isomorph_", rho_iso[i], T_iso[i], n, A_iso[i], "BoundedInversePower");
  }
}
