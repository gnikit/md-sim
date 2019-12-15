#pragma once
#include <algorithm>
#include <exception>
#include <numeric>
#include <tuple>
#include <vector>

/* A class containing useful methods that one can use */

class helper_functions {
 private:
  /* data */
 public:
  static std::vector<int> linspace(int a, int b, size_t N);

  static std::vector<double> linspace(double a, double b, size_t N);

  static std::tuple<double, double> linfit(const std::vector<double>& x,
                                           const std::vector<double>& y);
};