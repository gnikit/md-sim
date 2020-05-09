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
#include <stdlib.h>   /* abort */
#include <algorithm>  /* std::find */
#include <chrono>     /* CPU run-time */
#include <ctime>      /* std::chrono */
#include <functional> /* funciton pointers */
#include <iomanip>    /* setprecision */
#include <memory>     /* unique_ptr */
#include <numeric>    /* accumulate */
#include <random>     /* normal_dist */
#include <sstream>    /* stringstream */
#include <vector>     /* vectors */
// todo: change vector to std::valarray
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
  vector_3d<double> r;     /* Position Arrays */
  vector_3d<double> v;     /* Velocity Arrays */
  vector_3d<double> f;     /* Force arrays */
  vector_3d<double> Cv;    /* VAF arrays */
  vector_3d<double> MSD_r; /* used in MSD calculation */
  vector_3d<double> MSD;   /* MSD arrays */
  vector_3d<double> sf;    /* Structure factor k-arrays */

  /* Statistical quantity vectors, VAF, MSD, Energies and pressures */
  std::vector<double> Cr, msd, u_en, k_en, pc, pk, temperature, density, rdf,
      run_stats;

  /* Constructor variables */
  options_type options;

  /* Visualisation vectors, initialised in constructor */
  // todo: use vector of smart_pointers
  std::vector<std::vector<double>> *pos;

 private:
  double const PI = acos(-1.0);
  size_t step_idx = 0;
  /* Variables for storing inside the object the file ID */
  stat_file logger;

 public:
  MD();
  MD(options_type &input_options);
  MD(size_t step_number, std::vector<size_t> particles, std::string lattice);
  ~MD();

  /**
   * @brief
   * Executes the fluid simulation. It includes all statistical methods
   * defined in this class. It monitors the following quantities
   * RDF, MSD, VAF, average particle energies & pressures.
   * The produced files are named after the input parameters of the run.
   * The restart parameters of the run are also returned.
   */
  void simulation();

  /**
   * @brief
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
  void simulation(std::string const simulation_name, double const DENSITY,
                  double const TEMPERATURE, double POWER, double const A_CST,
                  std::string const pp_type);

  /**
   * @brief
   * Closes open file streams and resets sizes and values to 0
   * Use it when running multiple simulations and recycling the same
   * MD object.
   *
   * @param force_reset: Will close file streams even when compression is turned
   * on
   */
  void reset_values(bool force_reset = false);

  /**
   * @brief Fixes the random seed in the velocity generation distributions
   *
   * @param is_testing: boolean flag
   */
  void enable_testing(bool is_testing);

 protected:
  /*************************** INITIALISATION METHODS *************************/

  /**
   * @brief Initialises the:
   *  + Position Arrays
   *  + Velocity Arrays
   *  + Conserves/ Scales momentum
   *  + Temperature
   *  + Velocity Autocorrelation Function
   *
   * @param r: position vector of particles
   * @param v: velocity vectors of particles
   * @param TEMPERATURE: Thermostat target temperature
   * @return kinetic energy (normalised)
   */
  double initialise(vector_3d<double> &r, vector_3d<double> &v,
                    double const TEMPERATURE);

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
  void choose_lattice_formation(std::string const &lattice,
                                vector_3d<double> &r);

  /**
   * @brief
   * Generates velocities for based on the Maxwell Boltzmann distribution.
   * The MB dist in 3D is effectively the triple product of 3 Normal dist.
   *
   * The method populates the velocity vectors vx, vy, vz.
   * All the constants are assumed to be 1.
   *
   * @param v: velocity vectors of particles
   * @param TEMPERATURE: Temperature of the MB distribution
   */
  void mb_distribution(vector_3d<double> &v, double const TEMPERATURE);

  /**************************** ITERATIVE METHODS *****************************/

  /**
   * @brief Iterate the particles with one of the listed algorithms:
   *  - Explicit Verlet
   *  - Velocity Verlet
   *  - Runge Kutta 4th order
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param f: force vectors of particles
   * @param msd Flag set to true if calculating Mean Square Displacement.
   *            MSD requires a copy of the particles where BCs are not applied.
   * @return std::tuple<double, double, double> KE, U, PC
   */
  std::tuple<double, double, double> stepping_algorithm(vector_3d<double> &r,
                                                        vector_3d<double> &v,
                                                        vector_3d<double> &f,
                                                        bool msd);
  /**
   * @brief An iterative leap-frog Verlet Algorithm.
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param f: force vectors of particles
   * @return std::tuple<double, double, double> KE, U, PC
   */
  std::tuple<double, double, double> verlet(vector_3d<double> &r,
                                            vector_3d<double> &v,
                                            vector_3d<double> &f);

  /**
   * @brief Iterative Velocity verlet algorithm.
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param f: force vectors of particles
   * @return std::tuple<double, double, double> KE, U, PC
   */
  std::tuple<double, double, double> velocity_verlet(vector_3d<double> &r,
                                                     vector_3d<double> &v,
                                                     vector_3d<double> &f);

  /**
   * @brief A Runge Kutta 4th order iterative algorithm
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param f: force vectors of particles
   * @return std::tuple<double, double, double> KE, U, PC
   */
  std::tuple<double, double, double> runge_kutta4(vector_3d<double> &r,
                                                  vector_3d<double> &v,
                                                  vector_3d<double> &f);

  /**
   * @brief Calculates the forces interactions of a given pair potential
   * within the specified cutoff distance for a single step iteration
   *
   * @param r: position vectors of particles
   * @param f: force vectors of particles
   * @param potential: Type of the pair potential
   * @return std::tuple<double, double> Potential Energy, Configuration Pressure
   */
  std::tuple<double, double> calculate_forces(vector_3d<double> &x,
                                              vector_3d<double> &f,
                                              pair_potential_type &potential);

  /**************************** BOUNDARY CONDITIONS ***************************/

  /**
   * @brief Applies the boundary conditions:
   *  - Periodic
   *  - Hard Wall
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   * @param f: force vectors of particles
   */
  void apply_boundary_conditions(vector_3d<double> &r, vector_3d<double> &v,
                                 vector_3d<double> &f);

  /**
   * @brief Applies a periodic boundary condition, which simulates
   * an infinite volume of fluid.
   *
   * @param r: position vectors of particles
   */
  void periodic_boundary_conditions(vector_3d<double> &r);

  /**
   * @brief Applies a reflective wall boundary condition, which is
   * identical to having a hard wall boundary
   *
   * @param r: position vectors of particles
   * @param v: velocity vectors of particles
   */
  void reflective_boundary_conditions(vector_3d<double> const &r,
                                      vector_3d<double> &v);

  /********************** STATISTICAL QUANTITIES METHODS **********************/

  /**
   * @brief
   * Calculates the Velocity Autocorrelation Function for the fluid.
   * The method stores the values into the Cr vector and uses internally
   * the velocities of the particles.
   *
   * @param Cv: holds the particle velocities at t = 0
   * @param v: velocity vectors of particles
   */
  void velocity_autocorrelation_function(vector_3d<double> &Cv,
                                         vector_3d<double> &v);

  /**
   * @brief
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
   * @brief
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
  void mean_square_displacement(vector_3d<double> &MSD,
                                vector_3d<double> &MSD_r);

  /**
   * @brief
   * Calculates the structure factor, which is computed as a Fourier Transform.
   * The structure factor is more useful in FCC and BCC lattices
   *
   * @param r: position vectors of particles
   */
  void structure_factor(vector_3d<double> &r);

  /***************************** LOGGING METHODS ******************************/

  /**
   * @brief
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
  std::string set_simulation_params(double const &rho, double const &T,
                                    double const &power, double const &a,
                                    std::string const &pp_type);

  /**
   * @brief IO wrapper for saving into file the positions of the particles
   * at a given timestep. The file is a CSV with the following structure
   * for all particles.
   * x-pos, y-pos, z-pos, ...
   *
   * It relies on the flag options.io_options.visualise being set to true.
   */
  void save_visualisation_arrays(size_t dump_no);

  /**
   * @brief Manages all the IO operations, file streams, file creation etc.
   */
  void file_output(stat_file &logger);

  /**
   * @brief Get the statistics (min max mean l2norm rms) of the run quantities:
   *  + Mean Square Displacement
   *  + Velocity Autocorrelation Function
   *  + Structure Factor
   *  + Potential Energy
   *  + Configurational Pressure
   *
   * This routine is intended to be called at the end of the simulation run
   * so as the vector quantities have been correctly calculated.
   *
   * The quantities are not passed as arguments since the aim of this method is
   * to be called outside the MD class.
   *
   * @return std::vector<double> Vector containing all stats in the order listed
   */
  std::vector<double> calculate_run_stats();

  /**
   * @brief Get the run_stats vector that has been precalculated by
   *        calculate_run_stats
   *
   * @return std::vector<double> run_stats
   */
  std::vector<double> get_run_stats();

  /******************************* MISCELLANEOUS ******************************/

  /**
   * @brief Set the vector sizes. Resizes all vectors and containers
   */
  void set_vector_sizes();
};