
#pragma once
#include <tuple>

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

/*
 * A class responsible for implementing different
 * pair potentials in the MD class.
 *
 * Pair potentials supported:
 * BIP: Bounded Inverse Power
 * GCM: Gaussian Core Model
 */

class MD_tools {
 public:
  MD_tools();
  ~MD_tools();

  std::tuple<double, double> BIP_force(double &r, int n, double a);
  std::tuple<double, double> GCM_force(double &r);
  std::tuple<double, double> Exp_force(double &r, double m, double C);
};
