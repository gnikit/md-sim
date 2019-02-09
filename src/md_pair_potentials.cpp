#include "md_pair_potentials.h"

MD_tools::MD_tools() {}

MD_tools::~MD_tools() {}

std::tuple<double, double> MD_tools::get_force(double &r, double &power,
                                               double &a_cst,
                                               std::string &pp_name) {
  if (pp_name == "GCM") {
    auto [ff, u] = GCM_force(r);
    return std::make_tuple(ff, u);
  } else if (pp_name == "EXP") {
    auto [ff, u] = Exp_force(r, power, a_cst);
    return std::make_tuple(ff, u);
  } else {
    auto [ff, u] = BIP_force(r, power, a_cst);
    return std::make_tuple(ff, u);
  }
}

std::tuple<double, double> MD_tools::get_force(double &r) {
  auto [ff, u] = GCM_force(r);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::get_force(double &r, double &power,
                                               double &a_cst) {
  auto [ff, u] = BIP_force(r, power, a_cst);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::BIP_force(double &r, double n, double a) {
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
  double ff = (n)*r * pow(sqrt(r * r + a * a), ((-n - 2.0)));
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

std::tuple<double, double> MD_tools::Exp_force(double &r, double m, double C) {
  /*
   * Exponential pair potential, similar to the GCM.
   * If used use a smaller cut-off 1.5 ~ 2.0.
   * Also it is more convinient to use multipliers of e as the C cst
   */

  double ff = C * m * pow(r, (m - 1)) * exp(-(pow(r, m)));
  double u = C * exp(-(pow(r, m)));
  return std::make_tuple(ff, u);
}