#include "MD.h"
#include <chrono>   // CPU run-time
#include <ctime>    // std::chrono
#include <iomanip>  // setprecision
#include <random>   // normal_dist
#include <sstream>  // stringstream
#include "md_algorithms.h"

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

// FileIO has to be loaded after the math libraries
#include "FileIO.h"  // FileIO class

#define PARTICLES_PER_AXIS 10
// Increase the number of bins for a more accurate RDF
#define NHIST 500                // Number of histogram bins
#define RDF_WAIT 0               // Iterations after which RDF will be collected
#pragma warning(disable : 4996)  //_CRT_SECURE_NO_WARNINGS

MD::MD(std::string DIRECTORY, size_t run_number, bool COMPRESS_FLAG,
       size_t rdf_bins, size_t particles_per_axis, bool track_particles,
       size_t collect_rdf_after) {
  /*
   * This constructor allows for increased control over the internal
   * parameters of the fluid.
   */

  try {
    _dir = DIRECTORY;
    if (!fs::exists(_dir)) {
      // TODO: check the exception the fs::exists throws
      throw
          "input DIRECTORY in MD constructor does not exist. "
          "Use a valid directory for output files to be saved";
    }

  } catch (const char *msg) {
    std::cerr << "Error: " << msg << std::endl;
    exit(1);
  }

  STEPS = run_number;
  // Total number of particles N
  Nx = Ny = Nz = particles_per_axis;
  N = Nx * Ny * Nz;
  compression_flag = COMPRESS_FLAG;

  // Save all the positions for the fluid
  VISUALISE = track_particles;
  // Accuracy of RDF
  nhist = rdf_bins;

  // Ensuring the number of steps is greater than the rdf equilibration period
  try {
    // The number of iterations the data collection of RDF is postponed
    // in order to allow the fluid to lose its internal cubic lattice
    rdf_wait = collect_rdf_after;
    // Substraction of size_ts if negative results into logic errors
    // hence the use of an int temp;
    int temp = run_number - collect_rdf_after;
    if (temp < 0) {
      throw "collect_rdf_after is greater than the run_number";
    }
  } catch (const char *msg) {
    std::cerr << "Warning: " << msg << std::endl;
    std::cerr << "rdf_wait is set to 0" << std::endl;
    rdf_wait = 0;
  }
  // For efficiency, memory in the containers is reserved before use
  /* Positions */
  rx.reserve(N);
  ry.reserve(N);
  rz.reserve(N);
  /* Velocities */
  vx.reserve(N);
  vy.reserve(N);
  vz.reserve(N);
  /* RDF */
  gr.resize(nhist + 1, 0);  // gr with Index igr
  /* Forces/Acceleration */
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
  /* Observed Quantities */
  Cr.reserve(STEPS);    // Velocity Autocorrelation Function
  msd.reserve(STEPS);   // Mean Square Displacement
  u_en.reserve(STEPS);  // Average Potential Energy
  k_en.reserve(STEPS);  // Average Kinetic Energy
  pc.reserve(STEPS);    // Configurational Pressure
  pk.reserve(STEPS);    // Kinetic Pressure
  temperature.reserve(STEPS);
  /* Visualisation vectors on the heap*/
  pos_x = new std::vector<std::vector<double>>(STEPS);
  pos_y = new std::vector<std::vector<double>>(STEPS);
  pos_z = new std::vector<std::vector<double>>(STEPS);

  PI = acos(-1.0);
  _density_increment = 0.01;
  steps_quench = 10000;      // steps between each quenching

}

// Delegating constructors with reduced number of arguments
// https://en.wikipedia.org/wiki/C++11#Object_construction_improvement
// Constructor to use for density compression
MD::MD(std::string DIRECTORY, size_t run_number, bool COMPRESS_FLAG)
    : MD(DIRECTORY, run_number, COMPRESS_FLAG, NHIST, PARTICLES_PER_AXIS, false,
         RDF_WAIT) {}

// Constructor to use for the simplest cases
MD::MD(std::string DIRECTORY, size_t run_number)
    : MD(DIRECTORY, run_number, false, NHIST, PARTICLES_PER_AXIS, false,
         RDF_WAIT) {}

MD::~MD() {
  // Destroy the vectors allocated on the heap
  delete pos_x;
  delete pos_y;
  delete pos_z;
}

// Methods for MD Analysis
void MD::initialise(std::vector<double> &x, std::vector<double> &y,
                    std::vector<double> &z, std::vector<double> &vx,
                    std::vector<double> &vy, std::vector<double> &vz,
                    double TEMPERATURE) {
  /*
   *  Initialises the:
   *  + Position Arrays
   *  + Velocity Arrays
   *  + Conserves/ Scales momentum == 0
   *  + Temperature
   *  + Velocity Autocorrelation Function
   *
   *  @param &x, &y, &z: X, Y, Z vector points
   *  @param &vx, &vy, &vz: Vx, Vy, Vz vector points
   *  @param TEMPERATURE: Thermostat target temperature
   */

  // Initialise position matrix and velocity matrix from Cubic Centred Lattice
  if (compression_flag == false ||
      (compression_flag == true && Q_counter == 0)) {
    size_t n = 0;
    size_t i, j, k;
    for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) {
        for (k = 0; k < Nz; ++k) {
          x.push_back((i + 0.5) * scale);
          y.push_back((j + 0.5) * scale);
          z.push_back((k + 0.5) * scale);

          // Position vectors used for MSD
          rrx.push_back((i + 0.5) * scale);
          rry.push_back((j + 0.5) * scale);
          rrz.push_back((k + 0.5) * scale);

          ++n;
        }
      }
    }
    // Generates Maxwell-Boltzmann distribution
    mb_distribution(TEMPERATURE);
  }

  // scale of x, y, z
  double mean_vx = 0;
  double mean_vy = 0;
  double mean_vz = 0;

  size_t i;
  // Momentum conservation
  for (i = 0; i < N; ++i) {
    mean_vx += vx[i] / N;
    mean_vy += vy[i] / N;
    mean_vz += vz[i] / N;
  }

  size_t tempN = N;
  // Subtracting Av. velocities from each particle
  for (i = 0; i < tempN; ++i) {
    vx[i] = vx[i] - mean_vx;
    vy[i] = vy[i] - mean_vy;
    vz[i] = vz[i] - mean_vz;
  }
  // Temperature calculation, statistically
  KE = 0;
  for (i = 0; i < N; ++i) {
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  T = KE / (1.5 * N);
  scale_v = sqrt(TEMPERATURE / T);  // scaling factor

  // Velocity scaling
  for (i = 0; i < tempN; ++i) {
    vx[i] *= scale_v;
    vy[i] *= scale_v;
    vz[i] *= scale_v;
  }
  // MSD initialisation, storing first positions of particles
  MSDx = x;
  MSDy = y;
  MSDz = z;

  // VAF initialisation, storing first velocities of particles
  Cvx = vx;
  Cvy = vy;
  Cvz = vz;
}

void MD::mb_distribution(double TEMPERATURE) {
  /*
   * Generates velocities for based on the Maxwell Boltzmann distribution.
   * The MB dist in 3D is effectively the triple product of 3 Normal dist.
   *
   * The method populates the velocity vectors vx, vy, vz.
   * All the constants are assumed to be 1.
   *
   * @ param TEMPERATURE: Temperature of the MB distribution
   */

  double kb = 1.0;
  double m = 1.0;

  double var = sqrt(TEMPERATURE * kb / m);
  double mean = 0;

  // TODO: this needs a seed to randomise every time, FIX
  std::default_random_engine generator;
  std::normal_distribution<double> g_x(mean, var);
  std::normal_distribution<double> g_y(mean, var);
  std::normal_distribution<double> g_z(mean, var);

  for (size_t i = 0; i < N; ++i) {
    vx.push_back(g_x(generator));
    vy.push_back(g_y(generator));
    vz.push_back(g_z(generator));
  }
}

void MD::verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry,
                          std::vector<double> &rz, std::vector<double> &vx,
                          std::vector<double> &vy, std::vector<double> &vz,
                          std::vector<double> &rrx, std::vector<double> &rry,
                          std::vector<double> &rrz) {
  /*
   *  An iterative leap-frog Verlet Algorithm.
   *
   *  @params <r> (x,y,z): position vector of particles
   *  @params <v> (vx,vy,vz): velocity vector of particles
   *  @params <rr> (xx,yy,zz): position vectors for the fluid in the
   *                           square displacement arrays
   */

  size_t i;
  for (i = 0; i < N; ++i) {
    vx[i] = vx[i] * scale_v + fx[i] * dt;
    vy[i] = vy[i] * scale_v + fy[i] * dt;
    vz[i] = vz[i] * scale_v + fz[i] * dt;
    rx[i] = rx[i] + vx[i] * dt;
    ry[i] = ry[i] + vy[i] * dt;
    rz[i] = rz[i] + vz[i] * dt;

    rrx[i] = rrx[i] + vx[i] * dt;
    rry[i] = rry[i] + vy[i] * dt;
    rrz[i] = rrz[i] + vz[i] * dt;

    // Kinetic Energy Calculation
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    // Apply periodic boundary conditions to ensure particles remain
    // inside the box
    if (rx[i] > L) rx[i] = rx[i] - L;
    if (rx[i] < 0.0) rx[i] = rx[i] + L;
    if (ry[i] > L) ry[i] = ry[i] - L;
    if (ry[i] < 0.0) ry[i] = ry[i] + L;
    if (rz[i] > L) rz[i] = rz[i] - L;
    if (rz[i] < 0.0) rz[i] = rz[i] + L;
  }
}

void MD::velocity_autocorrelation_function(std::vector<double> &Cvx,
                                           std::vector<double> &Cvy,
                                           std::vector<double> &Cvz) {
  /*
   * Calculates the Velocity Autocorrelation Function for the fluid.
   * The method stores the values into the Cr vector and uses internally
   * the velocities of the particles.
   *
   * @param &Cvx, &Cvy, &Cvz: Reference to the initial velocity vectors.
   */

  double cr_temp = 0;  // resets the sum every time step
  double kb = 1.0;     // Boltzmann constant
  double m = 1.0;      // particle mass
  size_t i;
  /* The peak of the VAF is located at 3kb*T/m */
  double norm = 3 * _T0 / m;
  for (i = 0; i < N; i++) {
    cr_temp += (Cvx[i] * vx[i] + Cvy[i] * vy[i] + Cvz[i] * vz[i]);
  }
  Cr.push_back((cr_temp / N) / norm);
}

void MD::radial_distribution_function() {
  /*
   * Calculates the Radial Distribution Function for the fluid
   * based on the values that are stored in the gr vector
   * throughout the simulation.
   * The accuracy of the RDF is defined by the nhist variable.
   * Saves the unormalised and normalised values of the RDF in
   * the corresponding filestream.
   */

  double R = 0;
  double norm = 1;
  double cor_rho = _rho * (N - 1) / N;
  double dr = rg / nhist;
  size_t i;
  Hist << "# particles (N): " << N << " steps: " << STEPS << " rho: " << _rho
       << " bin (NHIST): " << nhist << " cut_off (rg): " << rg << " dr: " << dr
       << std::endl;
  Hist << "# Unormalised" << '\t' << "Normalised" << std::endl;
  for (i = 1; i < nhist; ++i) {  // Changed initial loop value from 0 -> 1
    R = rg * i / nhist;
    // Volume between 2 spheres, accounting for double counting
    // hence the 2/3*pi*((R+dr)**3 - R**3)
    // Accounting for the rdf_wait time steps
    norm = cor_rho * (2.0 / 3.0 * PI * N * (STEPS - rdf_wait) *
                      (pow((R + (dr / 2.0)), 3) - pow((R - (dr / 2.0)), 3)));

    Hist << gr[i] << '\t' << gr[i] / norm << std::endl;
  }
}

void MD::mean_square_displacement(std::vector<double> &MSDx,
                                  std::vector<double> &MSDy,
                                  std::vector<double> &MSDz) {
  /*
   * Performs the Mean Square Displacement calculation for the fluid.
   * The method stores the MSD of the step in the msd vector.
   *
   * The method uses the explicitly created position vectors
   * rrx, rry and rrz that are instantiated in the initialise method
   * and are a copy of the rx, ry, rz vectors without the application
   * of the periodic boundary conditions in the Verlet algorithm.
   *
   *
   * @param &MSDx, &MSDy, &MSDz: Reference to vectors containing the
   *                             initial positions of the particles
   */

  double msd_temp = 0;
  for (size_t i = 0; i < N; ++i) {
    msd_temp += (pow((rrx[i] - MSDx[i]), 2) + pow((rry[i] - MSDy[i]), 2) +
                 pow((rrz[i] - MSDz[i]), 2));
  }
  msd.push_back(msd_temp / N);
}

void MD::density_compression(int steps_quench, double TEMPERATURE,
                             double density_increment) {
  /*
   * Performs repeated compressions of the fluid by periodically
   * incrementing the density of the fluid.
   * As a consequence the box length, the scaling factor and the
   * position vectors are also scaled in order to conserve the number
   * of particles in the box.
   *
   * @param steps_quench: Number of steps over which a compression occurs
   * @param TEMPERATURE: Temperature of the fluid
   * @param density_increment: Amount with which the density is increased (or
   * decreased)
   */
  double old_box_length = L;

  // Density incrementation
  _rho += density_increment;

  // Re-using this piece of code from MD::Simulation
  scale = pow((N / _rho), (1.0 / 3.0)) / Nx;
  L = pow((N / _rho), 1.0 / 3.0);
  double box_length_ratio = L / old_box_length;

  // Rescalling the positional vectors
  for (size_t i = 0; i < N; ++i) {
    rx[i] *= box_length_ratio;
    ry[i] *= box_length_ratio;
    rz[i] *= box_length_ratio;
  }
  Vol = N / _rho;

  initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);
}

// MD Simulation
void MD::Simulation(double DENSITY, double TEMPERATURE, int POWER,
                    double A_CST) {
  /*
   *  Executes the fluid simulation. It includes all statistical methods
   *  defined in this class. It monitors the following quantities
   *  RDF, MSD, VAF, average particle energies & pressures.
   *
   *  @param DENSITY: Density rho of fluid.
   *  @param TEMPERATURE: Temperature of the thermostat.
   *  @param POWER: the power n that the pair potential will be using
   *                typical values are in the range of 6-18.
   *  @param A_CST: softening constant 'a' of the pair potential.
   *                When 'a' = 0, then fluid is pure MD, increasing
   *                'a' results into softening of the pair potential.
   */
  // Initialise scaling variables
  std::cout << "***************************\n"
               "** MD Simulation started **\n"
               "***************************\n"
            << std::endl;
  std::cout << "Fluid parameters: rho: " << DENSITY << " T: " << TEMPERATURE
            << " n: " << POWER << " a: " << A_CST << std::endl;
  _T0 = TEMPERATURE;
  _rho = DENSITY;
  dt /= sqrt(_T0);
  // Box length scaling
  scale = pow((N / _rho), (1.0 / 3.0)) / Nx;
  L = pow((N / _rho), 1.0 / 3.0);
  Vol = N / _rho;

  // cut_off redefinition
  // Large cut offs increase the runtime exponentially
  cut_off = 3.0;
  rg = cut_off;
  dr = rg / nhist;
  MD_tools potential;  // Pair potential object

  // Generating the filenames for the output
  data = file_naming("/Data", DENSITY, TEMPERATURE, POWER, A_CST);
  pos =
      file_naming("/Positions_Velocities", DENSITY, TEMPERATURE, POWER, A_CST);
  HIST = file_naming("/RDF", DENSITY, TEMPERATURE, POWER, A_CST);

  open_files();
  time_stamp(DATA, "# step \t rho \t U \t K \t Pc \t Pk \t MSD \t VAF");

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);

  for (_STEP_INDEX = 0; _STEP_INDEX < STEPS; ++_STEP_INDEX) {
    // Forces loop
    // Resetting forces
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    U = 0;  // seting Potential U to 0
    PC = 0;

    // TODO: pass as argument steps_quench
    if (compression_flag == true && _STEP_INDEX != 0 &&
        _STEP_INDEX % steps_quench == 0) {
      ++Q_counter;
      density_compression(steps_quench, TEMPERATURE, _density_increment);
      std::cout << "Compressing fluid, run: " << Q_counter << " rho: " << _rho
                << std::endl;
    }

    size_t i, j;
    for (i = 0; i < N - 1; ++i) {
      for (j = i + 1; j < N; ++j) {
        x = rx[i] - rx[j];  // Separation distance
        y = ry[i] - ry[j];  // between particles i and j
        z = rz[i] - rz[j];  // in Cartesian

        // Ensure the particles' separation distance is within the box
        if (x > (0.5 * L)) x = x - L;
        if (x < (-0.5 * L)) x = x + L;
        if (y > (0.5 * L)) y = y - L;
        if (y < (-0.5 * L)) y = y + L;
        if (z > (0.5 * L)) z = z - L;
        if (z < (-0.5 * L)) z = z + L;

        // Pair potential radius
        r = sqrt((x * x) + (y * y) + (z * z));

        // Force loop
        if (r < cut_off) {
          // BIP potential of the form: phi = 1/[(r**2 + a**2)**(n/2)]
          // Allows the user to choose different pair potentials
          auto [ff, temp_u] = potential.BIP_force(r, POWER, A_CST);

          // Average potential energy
          U += temp_u;

          // Configurational pressure
          PC += r * ff;

          // Canceling the ij and ji pairs
          // Taking the lower triangular matrix
          fx[i] += x * ff / r;
          fx[j] -= x * ff / r;
          fy[i] += y * ff / r;
          fy[j] -= y * ff / r;
          fz[i] += z * ff / r;
          fz[j] -= z * ff / r;

          // Radial Distribution
          // measured with a delay, since the system requires a few thousand
          // time-steps to reach equilibrium
          if (_STEP_INDEX > rdf_wait) {
            igr = round(nhist * r / rg);
            gr[igr] += 1;
          }
        }
      }
    }

    // Average Potential Energy per particle
    u_en.push_back(U / N);

    // Average Configurational Pressure Pc
    pc.push_back(PC / (3 * Vol));

    // Isothermal Calibration
    scale_v = sqrt(_T0 / T);  // using T & KE from prev timestep
    KE = 0;                   // resetting Kinetic Energy per iteration

    verlet_algorithm(rx, ry, rz, vx, vy, vz, rrx, rry, rrz);
    mean_square_displacement(MSDx, MSDy, MSDz);
    velocity_autocorrelation_function(Cvx, Cvy, Cvz);

    // Average Temperature
    T = KE / (1.5 * N);
    temperature.push_back(T);

    // Kinetic Pressure
    PK = _rho * T;
    pk.push_back(PK);

    // Average Kinetic Energy
    k_en.push_back(KE / N);

    // Density
    density.push_back(_rho);

    // Save positions for visualisation with Python
    if (VISUALISE) {
      // Reserve memory for the position vectors
      (*pos_x)[_STEP_INDEX].reserve(N);
      (*pos_y)[_STEP_INDEX].reserve(N);
      (*pos_z)[_STEP_INDEX].reserve(N);

      // Populate the vectors with the current positions
      (*pos_x)[_STEP_INDEX] = rx;
      (*pos_y)[_STEP_INDEX] = ry;
      (*pos_z)[_STEP_INDEX] = rz;
    }
  }
  // Simulation Ends HERE

  if (VISUALISE) {
    // Save particle positions to files
    FileIO f;
    // Write the arrays as jagged,(hence transposed), this creates rows=STEPS
    // and columns=PARTICLES
    f.Write2File<double>(
        *pos_x, file_naming("/x_data", DENSITY, TEMPERATURE, POWER, A_CST),
        "\t", true);
    f.Write2File<double>(
        *pos_y, file_naming("/y_data", DENSITY, TEMPERATURE, POWER, A_CST),
        "\t", true);
    f.Write2File<double>(
        *pos_z, file_naming("/z_data", DENSITY, TEMPERATURE, POWER, A_CST),
        "\t", true);
  }

  write_to_files();
  // Saving Last Position
  time_stamp(POS, "# X\tY\tZ\tVx\tVy\tVz\tFx\tFy\tFz");
  for (size_t el = 0; el < rx.size(); ++el) {
    POS << rx[el] << '\t' << ry[el] << '\t' << rz[el] << '\t' << vx[el] << '\t'
        << vy[el] << '\t' << vz[el] << '\t' << fx[el] << '\t' << fy[el] << '\t'
        << fz[el] << std::endl;
  }

  radial_distribution_function();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
      << "CPU run time = "
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() /
             60
      << " min "
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() %
             60
      << "s" << std::endl;
  std::cout << "******************************\n"
               "** MD Simulation terminated **\n"
               "******************************\n"
            << std::endl;
  // Close file streams and reset vectors
  reset_values();
}

// File Handling
std::string MD::file_naming(std::string prefix, double DENSITY,
                            double TEMPERATURE, int POWER, double A_cst) {
  /*
   * Generates a unique filename for the simulation results to be stored.
   * The method infers from the constructor the number of particles used
   * and the duration of the simulation (steps).
   *
   * @param prefix: File name preceeding the file id
   * @param DENSITY: Fluid density
   * @param TEMPERATURE: Fluid temperature
   * @param POWER: Pair potential power
   * @param A_cst: Softening parameter constant
   *
   * @return: string structured as follows
   *          INPUT_DIR/prefix_step_#_particles_#_rho_#_T_#_n_#_A_#.txt
   */

  // Individual streams handling double to string conversions
  std::stringstream A_stream, rho_stream, T_stream;

  rho_stream << std::fixed << std::setprecision(4) << DENSITY;    // 4 decimals
  T_stream << std::fixed << std::setprecision(4) << TEMPERATURE;  // 4 decimals
  A_stream << std::fixed << std::setprecision(5) << A_cst;        // 5 decimals

  _step_to_str = "_step_" + std::to_string(STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();
  _n_to_str = "_n_" + std::to_string(POWER);
  _A_to_str = "_A_" + A_stream.str();

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str + _T_to_str +
             _n_to_str + _A_to_str;

  // Explicit defitions
  _FILE_EXT = ".txt";

  // Path addition
  return _dir + prefix + _FILE_ID + _FILE_EXT;
}

void MD::open_files() {
  /*
   * Open/Create if file does not exist.
   * Overwrites existing data.
   */
  Hist.open(HIST, std::ios::out | std::ios::trunc);
  DATA.open(data, std::ios::out | std::ios::trunc);
  POS.open(pos, std::ios::out | std::ios::trunc);
}

void MD::write_to_files() {
  /*
   * Writes values of parameters to Data file.
   */
  for (size_t i = 0; i < STEPS; ++i) {
    DATA << (i + 1) << '\t' << density[i] << '\t' << temperature[i] << '\t'
         << u_en[i] << '\t' << k_en[i] << '\t' << pc[i] << '\t' << pk[i] << '\t'
         << msd[i] << '\t' << Cr[i] << std::endl;
  }
}

void MD::show_run(size_t step_size_show) {
  /*
   * Displays the system parameters every step_size_show of steps
   * Input the increment step.
   *
   * @param step_size_show: Number of steps, every which the parameters should
   * be displayed on screen
   */
  if (_STEP_INDEX == 0) {
    std::cout << "step:\tT:\tKE:\tU:\tU+K:\tPC:\tPK:\t(PK+PC):" << std::endl;
  }

  if (_STEP_INDEX % step_size_show == 0 || _STEP_INDEX == 1) {
    std::cout.precision(5);
    std::cout << _STEP_INDEX << "\t" << T << "\t" << KE << "\t" << U << "\t"
              << (U + KE) << "\t" << PC << "\t" << PK << "\t" << (PK + PC)
              << std::endl;
  }
}

void MD::reset_values() {
  /*
   * Closes open file streams and resets sizes and values to 0
   * Use it when running multiple simulations and recycling the same
   * MD object.
   */

  // Close streams
  Hist.close();
  DATA.close();
  POS.close();
  // Clear values, size, but reserve capacity
  rx.clear();
  ry.clear();
  rz.clear();
  rrx.clear();
  rry.clear();
  rrz.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  density.clear();
  temperature.clear();
  u_en.clear();
  k_en.clear();
  pc.clear();
  pk.clear();
  msd.clear();
  Cr.clear();
  gr.resize(nhist + 1, 0);  // gr with Index igr
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
}

void MD::time_stamp(std::ofstream &stream, std::string variables) {
  /*
   * Dates the file and allows the input of a header
   * Input a file stream to write and string of characters to display as
   * headers.
   *
   * @param &stream: Stream that should be timestaped
   * @param variables: A header that can be included underneath the timestamp
   */
  std::chrono::time_point<std::chrono::system_clock> instance;
  instance = std::chrono::system_clock::now();
  std::time_t date_time = std::chrono::system_clock::to_time_t(instance);
  stream << "# Created on: " << std::ctime(&date_time);
  stream << variables << std::endl;
}

std::string MD::convert_to_string(const double &x, const int &precision) {
  /*
   * Convert doubles to a string with a variable degree of precision.
   *
   * @param &x: Double number to be converted
   * @param &precision: Precision of the double when converted to string
   *
   * @return: string
   */
  std::ostringstream ss;
  ss.str(std::string());  // don't forget to empty the stream
  ss << std::fixed << std::setprecision(precision) << x;

  return ss.str();
}
