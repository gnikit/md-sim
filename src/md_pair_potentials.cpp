#include "md_pair_potentials.h"

std::tuple<double, double> MD_tools::BIP_force(double &r, double n, double a) {
  /* Force for particles separated a distance r */
  double ff = (n)*r * pow(sqrt(r * r + a * a), ((-n - 2.0)));
  double u = pow(sqrt(r * r + a * a), (-n));
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
                                             double a_cst) {
  auto [ff, u] = MD_tools::BIP_force(r, power, a_cst);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> GCM_pp::get_force(double &r, double m = NAN,
                                             double C = NAN) {
  auto [ff, u] = MD_tools::GCM_force(r);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> Exp_pp::get_force(double &r, double m, double C) {
  auto [ff, u] = MD_tools::Exp_force(r, m, C);
  return std::make_tuple(ff, u);
}

std::tuple<double, double> LJ_pp::get_force(double &r, double m = NAN,
                                            double C = NAN) {
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