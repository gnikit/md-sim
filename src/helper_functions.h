#pragma once
#include <algorithm>
#include <exception>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
/* Load Intel math lib if available */
#if defined(__INTEL_COMPILER)
#include <mathimf.h> /* Intel Math library */
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

/**
 * @brief Produces an equally spaced vector of N increments in the inclusive
 *        range of [a, b]
 *
 * @param a Initial range value
 * @param b Final range value (inclusive)
 * @param N Number of increments
 * @return std::vector<int> Equally spaced vector
 */
std::vector<int> linspace(int const& a, int const& b, size_t const& N);

/**
 * @brief Produces an equally spaced vector of N increments in the inclusive
 *        range of [a, b]
 *
 * @param a Initial range value
 * @param b Final range value (inclusive)
 * @param N Number of increments
 * @return std::vector<double> Equally spaced vector
 */
std::vector<double> linspace(double const& a, double const& b, size_t const& N);

/**
 * @brief Performs a linear regression fit on the x, y dataset provided
 *
 * @param x x-axis values
 * @param y y-axis values
 * @return std::tuple<double, double> slope, y intercept
 */
std::tuple<double, double> linfit(std::vector<double> const& x,
                                  std::vector<double> const& y);

/**
 * @brief Calculates the RMS of vector x
 *
 * @param x vector
 * @return double RMS
 */
double rms(std::vector<double> const& x);

/**
 * @brief Calculates the L2Norm of vector x
 *
 * @param x vector
 * @return double L2Norm
 */
double l2norm(std::vector<double> const& x);

/**
 * @brief Calculates the min, max, mean, l2norm and rms of vector x
 *
 * @param x vector
 * @return std::vector<double> vector of values
 */
std::vector<double> vector_stats(std::vector<double> const& x);

/**
 * @brief Repeats a string a supplied number of times
 *
 * @param str String to repeat
 * @param times Number of times to repeat
 * @return std::string Repeated string
 */
std::string repeat(std::string const& str, int times);

/**
 * @brief Pads a string with empty spaces
 *
 * @param str String to pad
 * @param pad Number of white spaces to add
 * @return std::string Padded string
 */
std::string pad_string(std::string const& str, size_t const& pad);
