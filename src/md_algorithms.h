
#pragma once
#include <tuple>

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"

#elif defined(__GNUC__)
#include <math.h>
#define COMPILER "G++"

#else
#include <math.h>
#define COMPILER "OTHER COMPILER"
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
};
