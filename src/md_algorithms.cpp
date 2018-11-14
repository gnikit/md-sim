#include "md_algorithms.h"

MD_tools::MD_tools() {
}

MD_tools::~MD_tools() {
}

std::tuple<double, double> MD_tools::BIP_force(double &r, int n, double a) {
  /*
   * Generates the force of a Bounded Inverse Power potential
   * and its potential energy.
   * 
   * @param &r: Reference to the separation distance between two particles
   * @param n: Pair potential strength
   * @param a: Softening parameter
   * 
   * @return std::tuple<double, double>: <Force, Potential energy>
   */
  // Force for particles separated a distance r
  double ff =
      (n)*r * pow(sqrt(r * r + a * a), ((-n - 2.0)));
  double u = pow(sqrt(r * r + a * a), (-n));
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::GCM_force(double &r) {
  /*
   * Generates the force of a Gaussian Core Model potential
   * and its potential energy.
   * 
   * @param &r: Reference to the separation distance between two particles
   * 
   * @return std::tuple<double, double>: <Force, Potential energy>
   */
  // Gaussian-force with sigma=1 and epsilon=1
  double ff = 2 * r * exp(-r * r);
  double u = exp(-r * r);
  return std::make_tuple(ff, u);
}
