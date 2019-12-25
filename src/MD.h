/**
 * Ioannis Nikiteas 13/7/2017
 *
 * BSc Dissertation:
 * Investigating the transition from Molecular Dynamics to
 * Smoothed Particle Hydrodynamics
 *
 * University: Royal Holloway University of London
 */
#pragma once
#include <algorithm>  /* std::find */
#include <chrono>     /* CPU run-time */
#include <ctime>      /* std::chrono */
#include <functional> /* funciton pointers */
#include <iomanip>    /* setprecision */
#include <numeric>    /* accumulate */
#include <random>     /* normal_dist */
#include <sstream>    /* stringstream */
#include <vector>     /* vectors */

#include "data_structures.h"
#include "helper_functions.h"
#include "md_pair_potentials.h"
#include "stat_file_logger.h"

/* Check for Compiler support */
// TODO: in future C++ versions, rm fs:: from global scope and mv in constructor
#if __cplusplus <= 201103L
#error This library requires at least C++17 compiler support
/* If C++ version C++2a or above use */
#elif __cplusplus >= 201709
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

/* Load Intel math lib if available */
#if defined(__INTEL_COMPILER)
#include <mathimf.h> /* Intel Math library */
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

class MD {
 protected:
  vector_3d r;     /* Position Arrays */
  vector_3d v;     /* Velocity Arrays */
  vector_3d f;     /* Force arrays */
  vector_3d Cv;    /* VAF arrays */
  vector_3d MSD_r; /* used in MSD calculation */
  vector_3d MSD;   /* MSD arrays */
  vector_3d sf;    /* Structure factor k-arrays */

  /* Statistical quantity vectors, VAF, MSD, Energies and pressures */
  std::vector<double> Cr, msd, u_en, k_en, pc, pk, temperature, density, rdf;

  /* Constructor variables */
  options_type options;

  /* Visualisation vectors, initialised in constructor */
  std::vector<std::vector<double>> *pos_x;
  std::vector<std::vector<double>> *pos_y;
  std::vector<std::vector<double>> *pos_z;

 private:
  double PI;
  /* Variables for storing inside the object the file ID */
  stat_file logger;

 public:
  MD();
  MD(options_type &input_options);
  MD(size_t step_number, std::vector<size_t> particles, std::string lattice);
  ~MD();

  /**
   * Executes the fluid simulation. It includes all statistical methods
   * defined in this class. It monitors the following quantities
   * RDF, MSD, VAF, average particle energies & pressures.
   * The produced files are named after the input parameters of the run.
   * The restart parameters of the run are also returned.
   */
  void simulation();

  /**
   * Identical to the simulation() but uses normal arguments
   * rather than the custom options data type. This results into
   * limited functionality.
   * @warning This method is not tested anymore
   *
   * @param DENSITY: Density rho of fluid.
   * @param TEMPERATURE: Temperature of the thermostat.
   * @param POWER: the power n that the pair potential will be using
   *               typical values are in the range of 6-18.
   * @param A_CST: softening constant 'a' of the pair potential.
   *               When 'a' = 0, then fluid is pure MD, increasing
   *               'a' results into softening of the pair potential.
   * @param pp_type: The type of the pair potential the simulation is
   *                 modelling. Options are "BoundedInversePower",
   *                 "GaussianCoreModel", "Exponential", "LennardJones"
   */
  void simulation(std::string simulation_name, double DENSITY,
                  double TEMPERATURE, double POWER, double A_CST,
                  std::string pp_type);

  /**
   * Closes open file streams and resets sizes and values to 0
   * Use it when running multiple simulations and recycling the same
   * MD object.
   *
   * @force_reset: Will close file streams even when compression is turned on
   */
  void reset_values(bool force_reset = false);

  /**
   * @brief Fixes the random seed in the velocity generation distributions
   *
   * @param is_testing: boolean flag
   */
  void enable_testing(bool is_testing);

 protected:
  /**
   * @brief Set the vector sizes. Resizes all vectors and containers
   *
   */
  void set_vector_sizes();

  /**
   * @brief Outputs the positions r of the particles, based on the input
   * lattice formation. The available lattice formations are
   *    + Simple Cubic
   *    + Face Centred Cubic
   *    + Body Centred Cubic
   *    + Random positions
   *
   * @param lattice: type of lattice
   * @param r: particle positions
   */
  void choose_lattice_formation(std::string &lattice, vector_3d &r);

  /**
   * Initialises the:
   * + Position Arrays
   * + Velocity Arrays
   * + Conserves/ Scales momentum
   * + Temperature
   * + Velocity Autocorrelation Function
   *
   * @param r: position vector of particles
   * @param v: velocity vectors of particles
   * @param TEMPERATURE: Thermostat target temperature
   * @return kinetic energy (normalised)
   */
  double initialise(vector_3d &r, vector_3d &v, double TEMPERATURE);

  /**
   * Generates velocities for based on the Maxwell Boltzmann distribution.
   * The MB dist in 3D is effectively the triple product of 3 Normal dist.
   *
   * The method populates the velocity vectors vx, vy, vz.
   * All the constants are assumed to be 1.
   *
   * @param v: velocity vectors of particles
   * @param TEMPERATURE: Temperature of the MB distribution
   */
  void mb_distribution(vector_3d &v, double TEMPERATURE);

  /**
   *  An iterative leap-frog Verlet Algorithm.
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param sample_msd: flag for calculating MSD, msd vectors are globally
   *                    defined and not passed as arguments
   * @return KE: kinetic energy (unormalised)
   */
  double verlet_algorithm(vector_3d &r, vector_3d &v, vector_3d &f, bool msd);

  double rk4_algorithm(vector_3d &r, vector_3d &v, vector_3d &f, bool msd);

  double stepping_algorithm(vector_3d &r, vector_3d &v, vector_3d &f, bool msd);

  /**
   * @brief Calculates the forces interactions of a given pair potential
   * within the specified cutoff distance for a single step iteration
   *
   * @param step_index: The number of the current iteration. Only used to
   *                    determine if the sampling of the RDF should occur, if
   *                    an equilibration period has been supplied.
   * @param potential: Type of the pair potential
   * @return std::tuple<double, double> Potential Energy, Configuration Pressure
   */
  std::tuple<double, double> calculate_forces(size_t &step_index,
                                              pair_potential_type potential);

  void apply_boundary_conditions();

  /**
   * Calculates the Velocity Autocorrelation Function for the fluid.
   * The method stores the values into the Cr vector and uses internally
   * the velocities of the particles.
   *
   * @param Cv: holds the particle velocities at t = 0
   * @param v: velocity vectors of particles
   */
  void velocity_autocorrelation_function(vector_3d &Cv, vector_3d &v);

  /**
   * Calculates the Radial Distribution Function for the fluid
   * based on the values that are stored in the gr vector
   * throughout the simulation.
   * The accuracy of the RDF is defined by the nhist variable.
   * Saves the unormalised and normalised values of the RDF in
   * the corresponding filestream.
   *
   * @param rho: density
   * @param cut_off: radius cutoff for which the simulation is evaluated
   * @param bins: number of bins, accuracy of the RDF histogram
   * @param particles: total number of particles
   */
  void radial_distribution_function(double &rho, double &cut_off, size_t &bins,
                                    size_t &particles, std::ofstream &fstream);

  /**
   * Performs the Mean Square Displacement calculation for the fluid.
   * The method stores the MSD of the step in the msd vector.
   *
   * The method uses the explicitly created position vectors
   * MSD_r that are instantiated in the initialise method
   * and are a copy of the r vectors without the application
   * of the periodic boundary conditions in the Verlet algorithm.
   *
   * @param MSD: vectors containing the particle positions at t = 0
   * @param MSD_r: position vectors of particles, with no boundaries
   */
  void mean_square_displacement(vector_3d &MSD, vector_3d &MSD_r);

  /**
   * Calculates the structure factor, which is computed as a Fourier Transform.
   * The structure factor is more useful in FCC and BCC lattices
   *
   * @param r: position vectors of particles
   */
  void structure_factor(vector_3d &r);

  /**
   * Sets the internal class variables for the fluid according to
   * the pair-potential selected.
   * Returns a string with all the simulation setup parameters
   *
   * @param rho: density
   * @param T: temperature
   * @param power: pair potential strength
   * @param a: softening parameter
   * @param pp_type: pair potential type
   * @return string containing simulation run parameters
   */ //todo: pass options type
  std::string set_simulation_params(double &rho, double &T, double &power,
                                    double &a, std::string &pp_type);

  /**
   * @brief IO wrapper for saving into file the positions of the particles
   * at each time step
   *
   */
  void save_visualisation_arrays();

  /**
   * @brief Manages all the IO operations, file streams, file creation etc.
   *
   */
  void file_output();
};