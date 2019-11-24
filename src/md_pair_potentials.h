
#pragma once
#include <iostream>
#include <map>
#include <tuple>

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

// *NOTE: currently the pair potentials are limited to a 3 argument parameters.
using pair_potential_type = std::tuple<double, double> (*)(double &, double,
                                                           double);

/*
 * A class responsible for implementing different
 * pair potentials in the MD class.
 *
 * Pair potentials supported:
 * BIP: Bounded Inverse Power
 * GCM: Gaussian Core Model
 * Exp: Exponential potential
 * LJ:  Lennard-Jones potential
 */

class MD_tools {
 public:
  static std::tuple<double, double> BIP_force(double &r, double n, double a);
  static std::tuple<double, double> GCM_force(double &r);
  static std::tuple<double, double> Exp_force(double &r, double m, double C);
  static std::tuple<double, double> LJ_force(double &r);
};

/* **********************  Pair potential wrappers  ************************ */

class BIP_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double n, double a);
};

class GCM_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m, double C);
};

class Exp_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m, double C);
};

class LJ_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m, double C);
};

pair_potential_type get_force_func(std::string pp_type);