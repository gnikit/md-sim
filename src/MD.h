/*
 * Ioannis Nikiteas 13/7/2017
 *
 * BSc Dissertation:
 * Investigating the transition from Molecular Dynamics to
 * Smoothed Particle Hydrodynamics
 *
 * University: Royal Holloway University of London
 *
 * A program meant to simulate a MD fluid with an only
 * repulsive BIP pair-potential. Increasing the parameter A
 * creates a coarse-graining effect for the system allowing it
 * to transition to SPH
 */
#pragma once
#include <algorithm>   // std::find
#include <chrono>      // CPU run-time
#include <ctime>       // std::chrono
#include <functional>  // funciton pointers
#include <iomanip>     // setprecision
#include <numeric>     // accumulate
#include <random>      // normal_dist
#include <sstream>     // stringstream
#include <vector>      // vectors

#include "md_pair_potentials.h"
#include "stat_file_logger.h"

// Check for Compiler support
// TODO: in future C++ versions, rm fs:: from global scope and mv in constructor
#if __cplusplus <= 201103L
#error This library requires at least C++17 compiler support
// If C++ version C++2a or above use
#elif __cplusplus >= 201709
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

// Load Intel math lib if available
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

class MD {
 protected:
  std::vector<double> rx, ry, rz;        // Position Arrays
  std::vector<double> vx, vy, vz;        // Velocity Arrays
  std::vector<double> fx, fy, fz;        // Force arrays
  std::vector<double> Cvx, Cvy, Cvz;     // VAF arrays
  std::vector<double> rrx, rry, rrz;     // used in MSD calculation
  std::vector<double> MSDx, MSDy, MSDz;  // MSD arrays
  std::vector<double> sfx, sfy, sfz;     // Structure factor k-arrays

  // Statistical quantity vectors, VAF, MSD, Energies and pressures
  std::vector<double> Cr, msd, u_en, k_en, pc, pk, temperature, density;

  /* Constructor variables */
  size_t __Nx, __Ny, __Nz, __N;  // Particles in x, y, z and total
  size_t __step_idx, __steps;    // step index, maximum steps
  size_t __nhist, __rdf_wait;    // histogram bin number, equilibrium period
  std::string __lattice;         // Parameter defining the initial lattice

  double __dt;       // time step
  double __KE;       // Kinetic Energy
  double __T;        // Temperature
  double __L;        // Length of the box after scaling
  double __cut_off;  // simulation runs only within cutoff
  double __var;      // variance random positions generator
  double scale_v;    // velocity scaling

  /* Visualisation vectors, initialised in constructor */
  std::vector<std::vector<double>> *pos_x;
  std::vector<std::vector<double>> *pos_y;
  std::vector<std::vector<double>> *pos_z;

  /* Compression variables */
  bool __compress;
  size_t c_counter = 0;  // counts the number compression that have occurred

  /* Radial Distribution variables */
  int igr;                 // Index of Hist
  double __rg;             // cut off radius
  double __dr;             // bin increment in terms of radius units
  std::vector<double> gr;  // RDF vector container

 private:
  double PI;
  /* Variables for storing inside the object the file ID */
  stat_file logger;
  std::string __dir;
  std::string __simulation_name;

  /* Variables storing the simulation setup variables */
  double __rho;           // density
  double __T0;            // target/ Thermostat temperature
  double __power;         // pair potential intensity
  double __a_cst;         // generic softening parameter
  std::string __pp_type;  // pair potential type

  /* Visualisation flag */
  bool __visualise;
  bool fixed_seed;

 public:
  /**
   * Contructor to initialise the fluid parameters
   *
   * @param step_number: number of iterations
   * @param particles: number of particles for x, y and z axes
   * @param lattice: type of lattice for the initial position of the particles
   *                 + Simple Cubic
   *                 + Face Centred Cubic
   *                 + Body Centred Cubic
   */
  MD(size_t step_number, std::vector<size_t> particles, std::string lattice);

  ~MD();

  /**
   * Load secondary options for a simulation run that allow for
   * better fine tuning of the model.
   *
   * @param out_directory: output directory of log files.
   * @param is_compressing: sets the commpression flag.
   * @param track_particles: saves the x, y, z positions of particles to
   *                         visualise them with python.
   * @param rdf_bins: the accuracy of the RDF.
   * @param collect_rdf_after: delay the collection of the RDF data to ensure
   *                           melting of the fluid.
   */
  void load_options(std::string out_directory, bool track_particles,
                    size_t rdf_bins, size_t collect_rdf_after);

  /**
   * Executes the fluid simulation. It includes all statistical methods
   * defined in this class. It monitors the following quantities
   * RDF, MSD, VAF, average particle energies & pressures.
   * The produced files are named after the input parameters of the run.
   * The restart parameters of the run are also returned.
   *
   * @param DENSITY: Density rho of fluid.
   * @param TEMPERATURE: Temperature of the thermostat.
   * @param POWER: the power n that the pair potential will be using
   *               typical values are in the range of 6-18.
   * @param A_CST: softening constant 'a' of the pair potential.
   *               When 'a' = 0, then fluid is pure MD, increasing
   *               'a' results into softening of the pair potential.
   * @param pp_type: The type of the pair potential the simulation is
   *                 modelling. Options are "BIP", "GCM", "EXP", "LJ"
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

 protected:
  /**
   * Initialises the:
   * + Position Arrays
   * + Velocity Arrays
   * + Conserves/ Scales momentum == 0
   * + Temperature
   * + Velocity Autocorrelation Function
   *
   * @param x, y, z: position vector of particles
   * @param vx, vy, vz: velocity vectors of particles
   * @param TEMPERATURE: Thermostat target temperature
   * @return kinetic energy (normalised)
   */
  double initialise(std::vector<double> &x, std::vector<double> &y,
                    std::vector<double> &z, std::vector<double> &vx,
                    std::vector<double> &vy, std::vector<double> &vz,
                    double TEMPERATURE);

  /**
   * Generates velocities for based on the Maxwell Boltzmann distribution.
   * The MB dist in 3D is effectively the triple product of 3 Normal dist.
   *
   * The method populates the velocity vectors vx, vy, vz.
   * All the constants are assumed to be 1.
   *
   * @param vx, vy, vz: velocity vectors of particles
   * @param TEMPERATURE: Temperature of the MB distribution
   */
  void mb_distribution(std::vector<double> &vx, std::vector<double> &vy,
                       std::vector<double> &vz, double TEMPERATURE);

  /**
   *  An iterative leap-frog Verlet Algorithm.
   *
   * @param x, y, z: position vectors of particles
   * @param vx, vy, vz: velocity vectors of particles
   * @param sample_msd: flag for calculating MSD, msd vectors are globally
   *                    defined and not passed as arguments
   * @return KE: kinetic energy (unormalised)
   */
  double verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry,
                          std::vector<double> &rz, std::vector<double> &vx,
                          std::vector<double> &vy, std::vector<double> &vz,
                          bool sample_msd);

  /**
   * Calculates the Velocity Autocorrelation Function for the fluid.
   * The method stores the values into the Cr vector and uses internally
   * the velocities of the particles.
   *
   * @param Cvx, Cvy, Cvz: holds the particle velocities at t = 0
   * @param vx, vy, vz: velocity vectors of particles
   */
  void velocity_autocorrelation_function(std::vector<double> &Cvx,
                                         std::vector<double> &Cvy,
                                         std::vector<double> &Cvz,
                                         std::vector<double> &vx,
                                         std::vector<double> &vy,
                                         std::vector<double> &vz);

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
                                    size_t &particles);

  /**
   * Performs the Mean Square Displacement calculation for the fluid.
   * The method stores the MSD of the step in the msd vector.
   *
   * The method uses the explicitly created position vectors
   * rrx, rry and rrz that are instantiated in the initialise method
   * and are a copy of the rx, ry, rz vectors without the application
   * of the periodic boundary conditions in the Verlet algorithm.
   *
   * @param MSDx, MSDy, MSDz: vectors containing the particle positions at t = 0
   * @param rrx, rry, rrz: position vectors of particles, with no boundaries
   */
  void mean_square_displacement(std::vector<double> &MSDx,
                                std::vector<double> &MSDy,
                                std::vector<double> &MSDz,
                                std::vector<double> &rrx,
                                std::vector<double> &rry,
                                std::vector<double> &rrz);

  /**
   * Calculates the structure factor, which is computed as a Fourier Transform.
   * The structure factor is more useful in FCC and BCC lattices
   *
   * @param rx, ry, rz: position vectors of particles
   */
  void structure_factor(std::vector<double> &rx, std::vector<double> &ry,
                        std::vector<double> &rz);

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
   */
  std::string set_simulation_params(double &rho, double &T, double &power,
                                    double &a, std::string &pp_type);
  /* Helper Functions */
 public:
  std::string get_dir();

  std::string get_simulation_name();

  std::string get_lattice_structure();

  bool get_visualisation_flag();

  void set_visualisation_flag(bool is_visualising);

  size_t get_particle_number();

  size_t get_rdf_accuracy();

  void enable_testing(bool is_testing);

  size_t get_rdf_collect_after();

  void set_rdf_collect_after(size_t rdf_collect_after);

  /**
   * Sets the variance for 3 normal distributions that will in turn
   * generate random particle positions in space.
   * To be called before the md_distribution function.
   * The md_distribution takes the square root of the Temperature hence the
   * square of the variance is being stored interlay.
   *
   * @param var: variance of normal distributions
   */
  void set_random_position_variance(double var);
};