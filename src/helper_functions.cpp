#include "helper_functions.h"

std::vector<int> helper_functions::linspace(int a, int b, size_t N) {
  /*
   * Produces an equally spaced vector of N increments
   * in the inclusive range of [a, b]
   */
  int h = (b - a) / (N - 1);
  std::vector<int> xs(N);
  std::vector<int>::iterator x;
  int val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
    *x = val;
  }
  return xs;
}

std::vector<double> helper_functions::linspace(double a, double b, size_t N) {
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
