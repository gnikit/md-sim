#include "md_pair_potentials.h"

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
  /* Force for particles separated a distance r */
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
  /* Gaussian-force with sigma=1 and epsilon=1 */
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

std::tuple<double, double> MD_tools::LJ_force(double &r) {
  /*
   * Lennard-Jones pair potential.
   */

  double ff = 4 * (12 * pow(r, -13) + 6 * pow(r, -7));
  double u = 4 * (pow(r, -12) - pow(r, -6));
  return std::make_tuple(ff, u);
}

/* **********************  Pair potential classes  ************************ */

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

/*
 * NOTE: this map needs to be updated with every new pair potential
 * that is added in the system.
 * In addition to the updating the map, a method has to be added in the MD_tools
 * and a wrapper class has to be created.
 */

std::map<std::string, pair_potential_type> get_force_funcs = {
    {"GCM", &GCM_pp::get_force}, {"GaussianCoreModel", &GCM_pp::get_force},
    {"EXP", &Exp_pp::get_force}, {"Exponential", &Exp_pp::get_force},
    {"BIP", &BIP_pp::get_force}, {"BoundedInversePower", &BIP_pp::get_force},
    {"LJ", &LJ_pp::get_force},   {"LennardJones", &LJ_pp::get_force}};

/* Global hash table function, returning a different type of object per key */
pair_potential_type get_force_func(std::string pp_type) {
  pair_potential_type pp = get_force_funcs[pp_type];
  if (!pp) pp = get_force_funcs["BIP"];  /* Default case */
  return pp;
}