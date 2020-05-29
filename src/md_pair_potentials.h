/**
 * @file md_pair_potentials.h
 * @author your name (you@domain.com)
 * @brief Contais the different pair potentials availble and their derivatives
 * for the calculation of the forces.
 * @version 0.1
 * @date 2020-05-10
 *
 * @copyright Copyright (c) 2020
 *
 */
#pragma once
#include <iostream>
#include <map>
#include <tuple>

/* Load Intel math lib if available */
#if defined(__INTEL_COMPILER)
#include <mathimf.h> /* Intel Math library */
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

/**
 * @brief Currently limited to pair potentials using at max 4 input parameters
 *
 */
using pair_potential_type = std::tuple<double, double> (*)(double &, double,
                                                           double, double);

/**
 * @brief
 * A class responsible for implementing different
 * pair potentials in the MD class.
 *
 * Pair potentials supported:
 * BIP or BoundedInversePower Bounded Inverse Power
 * GCM or GaussianCoreModel: Gaussian Core Model
 * EXP or Exponential: Exponential potential
 * LJ or LennardJones:  Lennard-Jones potential
 */
class MD_tools {
 public:
  /**
   * @brief
   * Generates the force of a Bounded Inverse Power potential
   * and its potential energy.
   *
   * @f[
   *  \phi(r) = \frac{1}{r^{q} + a^{q}}^{n/q}
   * @f]
   *
   * @f[
   *  F(r) = - \frac{\partial \phi}{\partial r} =
   *       n r^{q-1} (r^{q} + a^{q})^{-n/q - 1}
   * @f]
   *
   * @param r Separation distance between two particles
   * @param n Pair potential strength
   * @param a Softening parameter
   * @param q Moderating pair potential strength (generalised form)
   * @return std::tuple<double, double>: <Force, Potential energy>
   */
  static std::tuple<double, double> BIP_force(double &r, double n, double a,
                                              double q);

  /**
   * @brief
   * Generates the force of a Gaussian Core Model potential
   * and its potential energy.
   *
   * @param r Separation distance between two particles
   * @return std::tuple<double, double> <Force, Potential energy>
   */
  static std::tuple<double, double> GCM_force(double &r);

  /**
   * @brief Exponential pair potential, similar to the GCM.
   * If used use a smaller cut-off 1.5 ~ 2.0.
   * Also it is more convinient to use multipliers of e as the C constant
   *
   * @param r Separation distance between two particles
   * @param m Exponential parameter
   * @param C Scaling parameters
   * @return std::tuple<double, double>
   */
  static std::tuple<double, double> Exp_force(double &r, double m, double C);

  /**
   * @brief Lennard-Jones force calculation.
   *
   * @param r
   * @return std::tuple<double, double>
   */
  static std::tuple<double, double> LJ_force(double &r);
};

/*************************  Pair potential classes  ***************************/

class BIP_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double n, double a,
                                              double q = 2);
};

class GCM_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m = NAN,
                                              double C = NAN, double q = NAN);
};

class Exp_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m, double C,
                                              double q = NAN);
};

class LJ_pp {
 public:
  static std::tuple<double, double> get_force(double &r, double m = NAN,
                                              double C = NAN, double q = NAN);
};

/**************  Query Function for Pair Potential Hash Table *****************/

/**
 * @brief Get the force function for a given pair potential. It queries a local
 * (available only in md_pair_potentials.cpp) hash table to fetch the results.
 *
 * @param pp_type a string containing the pair potential name
 * @return pair_potential_type a tuple with pair potential quantities
 */
pair_potential_type get_force_func(std::string pp_type);