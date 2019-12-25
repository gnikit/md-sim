#pragma once
#include <algorithm>
#include <exception>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

/* A class containing useful methods that one can use */

class helper_functions {
 private:
  /* data */
 public:
  static std::vector<int> linspace(int const& a, int const& b, size_t const& N);

  static std::vector<double> linspace(double const& a, double const& b,
                                      size_t const& N);

  static std::tuple<double, double> linfit(std::vector<double> const& x,
                                           std::vector<double> const& y);

  static std::string repeat(std::string const& str, int times);

  static std::string pad_string(std::string const& str, size_t const& pad);
};