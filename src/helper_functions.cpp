#include "helper_functions.h"

std::vector<int> helper_functions::linspace(int const& a, int const& b,
                                            size_t const& N) {
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

std::vector<double> helper_functions::linspace(double const& a, double const& b,
                                               size_t const& N) {
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

std::tuple<double, double> helper_functions::linfit(
    std::vector<double> const& x, std::vector<double> const& y) {
  /*
   * Calculates a simple linear regression fit
   * and returns the slope (a) and y-intercept (b) of the fit.
   */

  /* Test the supplied vectors are of equal length */
  if (x.size() != y.size()) throw "X and Y are unequal sizes";

  const auto n = x.size();
  const auto s_x = std::accumulate(x.begin(), x.end(), 0.0);
  const auto s_y = std::accumulate(y.begin(), y.end(), 0.0);
  const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

  const auto denominator = (n * s_xx - s_x * s_x);
  if (denominator == 0) throw "Denominator is 0";

  /* Slope of the linear fit */
  const auto a = (n * s_xy - s_x * s_y) / denominator;
  /* Y-intercept of the linear fit */
  const auto b = (s_y * s_xx - s_x * s_xy) / denominator;

  return std::make_tuple(a, b);
}

std::string helper_functions::repeat(std::string const& str, int times) {
  std::string result;
  result.reserve(times * str.length());  // avoid repeated reallocation
  for (int a = 0; a < times; a++) result += str;
  return result;
}

std::string helper_functions::pad_string(std::string const& str,
                                         size_t const& pad) {
  size_t pad_size = 0;
  if (pad > str.size())
    pad_size = pad - str.size();
  else
    std::cerr << "Error: Attempting to pad a string with a smaller pad size"
                 " the actual string length"
              << std::endl;

  return str + repeat(" ", pad_size);
}