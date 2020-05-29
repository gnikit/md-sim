#include "md_pair_potentials.h"

std::tuple<double, double> MD_tools::BIP_force(double &r, double n, double a,
                                               double q) {
  /* Force for particles separated a distance r */

  double ar, ff, u;
  // The sqrt function is blazing fast there is no reason to not use it if q=2
  if (abs(q - 2.0) <= std::numeric_limits<double>::epsilon()) {
    ar = sqrt(r * r + a * a);
    ff = n * r * pow(ar, (-n - 2.0));
    u = pow(ar, -n);

  } else {
    ar = pow(r, q) + pow(a, q);
    ff = n * pow(r, q - 1.0) * pow(ar, (-n / q - 1.0));
    u = pow(ar, (-n / q));
  }
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::GCM_force(double &r) {
  /* Gaussian-force with sigma=1 and epsilon=1 */
  double ff = 2 * r * exp(-r * r);
  double u = exp(-r * r);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::Exp_force(double &r, double m, double C) {
  double ff = C * m * pow(r, (m - 1)) * exp(-(pow(r, m)));
  double u = C * exp(-(pow(r, m)));
  return std::make_tuple(ff, u);
}

std::tuple<double, double> MD_tools::LJ_force(double &r) {
  double ff = 4 * (12 * pow(r, -13) + 6 * pow(r, -7));
  double u = 4 * (pow(r, -12) - pow(r, -6));
  return std::make_tuple(ff, u);
}

/*************************  Pair potential classes  ***************************/

std::tuple<double, double> BIP_pp::get_force(double &r, double power,
                                             double a_cst, double q) {
  auto [ff, u] = MD_tools::BIP_force(r, power, a_cst, q);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> GCM_pp::get_force(double &r, double m, double C,
                                             double q) {
  auto [ff, u] = MD_tools::GCM_force(r);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> Exp_pp::get_force(double &r, double m, double C,
                                             double q) {
  auto [ff, u] = MD_tools::Exp_force(r, m, C);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> LJ_pp::get_force(double &r, double m, double C,
                                            double q) {
  auto [ff, u] = MD_tools::LJ_force(r);
  return std::make_tuple(ff, u);
}

/**************  Query Function for Pair Potential Hash Table *****************/

/**
 * @brief this map needs to be updated with every new pair potential
 * that is added in the system. In addition to the updating the map, a method
 * has to be added in the MD_tools and a wrapper class has to be created.
 *
 */
std::map<std::string, pair_potential_type> get_force_funcs = {
    {"GCM", &GCM_pp::get_force}, {"GaussianCoreModel", &GCM_pp::get_force},
    {"EXP", &Exp_pp::get_force}, {"Exponential", &Exp_pp::get_force},
    {"BIP", &BIP_pp::get_force}, {"BoundedInversePower", &BIP_pp::get_force},
    {"LJ", &LJ_pp::get_force},   {"LennardJones", &LJ_pp::get_force}};

pair_potential_type get_force_func(std::string pp_type) {
  pair_potential_type pp = get_force_funcs[pp_type];
  if (!pp) pp = get_force_funcs["BIP"]; /* Default case */
  return pp;
}