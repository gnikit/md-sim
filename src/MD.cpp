#include "MD.h"
#include <chrono>      // CPU run-time
#include <ctime>       // std::chrono
#include <functional>  // funciton pointers
#include <iomanip>     // setprecision
#include <numeric>     // accumulate
#include <random>      // normal_dist
#include <sstream>     // stringstream
#include "md_pair_potentials.h"

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

MD::MD(std::string &DIRECTORY, size_t run_number, bool COMPRESS_FLAG,
       size_t rdf_bins, size_t particles_per_axis, std::string LATTICE,
       bool track_particles, size_t collect_rdf_after) {
  /*
   * This constructor allows for increased control over the internal
   * parameters of the fluid.
   *
   * If the compress parameter is true, then the number of steps per simulation
   * automatically becomes the number of steps per compression.
   * The true duration of the simulation is the controlled by the initial and
   * final DENISITIES passed on the get_phases method and with their increment.
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
  if (LATTICE == "FCC") {
    lattice = 1;
    Nx = Ny = Nz = particles_per_axis;
    N = Nx * Ny * Nz * 4;
  } else if (LATTICE == "BCC")
    lattice = 2;
  else {
    lattice = 0;
    Nx = Ny = Nz = particles_per_axis;
    N = Nx * Ny * Nz;
  }

  // Initialise the log output variable
  logger.N = N;
  logger.STEPS = STEPS;
  logger._dir = _dir;

  // If compress is true, then STEPS = steps_per_compression
  compress = COMPRESS_FLAG;

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
  /* Structure factor k-arrays */
  sfx.reserve(N);
  sfy.reserve(N);
  sfz.reserve(N);
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
}

// Delegating constructors with reduced number of arguments
// https://en.wikipedia.org/wiki/C++11#Object_construction_improvement
// Constructor to use for density compress
MD::MD(std::string &DIRECTORY, size_t run_number, bool COMPRESS_FLAG)
    : MD(DIRECTORY, run_number, COMPRESS_FLAG, NHIST, PARTICLES_PER_AXIS, "SC",
         false, RDF_WAIT) {}

// Constructor to use for the simplest cases
MD::MD(std::string &DIRECTORY, size_t run_number)
    : MD(DIRECTORY, run_number, false, NHIST, PARTICLES_PER_AXIS, "SC", false,
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
  if (!compress || (compress && c_counter == 0)) {
    switch (lattice) {
      // FCC lattice
      case 1: {
        // Coordinates for the FCC lattice
        double x_c[4] = {0.25, 0.75, 0.75, 0.25};
        double y_c[4] = {0.25, 0.75, 0.25, 0.75};
        double z_c[4] = {0.25, 0.25, 0.75, 0.75};

        // Loop over the the corner coordinates of the FCC and then x, y, z
        for (size_t c = 0; c < 4; ++c) {
          for (size_t i = 0; i < Nx; ++i) {
            for (size_t j = 0; j < Ny; ++j) {
              for (size_t k = 0; k < Nz; ++k) {
                x.push_back((i + x_c[c]) / double(Nx) * scale);
                y.push_back((j + y_c[c]) / double(Ny) * scale);
                z.push_back((k + z_c[c]) / double(Nz) * scale);
              }
            }
          }
        }
        break;
      }

      // Simple Cubic lattice
      default: {
        for (size_t i = 0; i < Nx; ++i) {
          for (size_t j = 0; j < Ny; ++j) {
            for (size_t k = 0; k < Nz; ++k) {
              x.push_back((i + 0.5) * scale);
              y.push_back((j + 0.5) * scale);
              z.push_back((k + 0.5) * scale);
            }
          }
        }
        break;
      }
    }
    // Generates Maxwell-Boltzmann distribution
    mb_distribution(TEMPERATURE);
  }

  // Calculate the average velocities
  double mean_vx = std::accumulate(vx.begin(), vx.end(), 0.0) / N;
  double mean_vy = std::accumulate(vy.begin(), vy.end(), 0.0) / N;
  double mean_vz = std::accumulate(vz.begin(), vz.end(), 0.0) / N;
  // Conserve the momentum of the fluid by subsracting the average velocities
  // using a lambda expression
  std::for_each(vx.begin(), vx.end(), [mean_vx](double &d) { d -= mean_vx; });
  std::for_each(vy.begin(), vy.end(), [mean_vy](double &d) { d -= mean_vy; });
  std::for_each(vz.begin(), vz.end(), [mean_vz](double &d) { d -= mean_vz; });

  // ! Check the values of mean_vx, mean_vy, mean_vz after a compression
  size_t tempN = N;
  size_t i;
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
  // ! Check the value of KE after a compression
  T = KE / (1.5 * N);
  scale_v = sqrt(TEMPERATURE / T);  // scaling factor

  // Velocity scaling
  for (i = 0; i < tempN; ++i) {
    vx[i] *= scale_v;
    vy[i] *= scale_v;
    vz[i] *= scale_v;
  }
  // A copy of the r vectors where the BC will not be applied
  rrx = x;
  rry = y;
  rrz = z;
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

  // Use current time as seed for random generator
  std::srand(std::time(nullptr));
  int random_variable = std::rand();

  std::default_random_engine generator;
  generator.seed(random_variable);

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
  logger.RDF << "# particles (N): " << N << " steps: " << STEPS
             << " rho: " << _rho << " bin (NHIST): " << nhist
             << " cut_off (rg): " << rg << " dr: " << dr << std::endl;
  logger.RDF << "# Unormalised" << '\t' << "Normalised" << std::endl;
  for (i = 1; i < nhist; ++i) {  // Changed initial loop value from 0 -> 1
    R = rg * i / nhist;
    // Volume between 2 spheres, accounting for double counting
    // hence the 2/3*pi*((R+dr)**3 - R**3)
    // Accounting for the rdf_wait time steps
    norm = cor_rho * (2.0 / 3.0 * PI * N * (STEPS - rdf_wait) *
                      (pow((R + (dr / 2.0)), 3) - pow((R - (dr / 2.0)), 3)));

    logger.RDF << gr[i] << '\t' << gr[i] / norm << std::endl;
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
void MD::structure_factor(std::vector<double> &rx, std::vector<double> &ry,
                          std::vector<double> &rz) {
  /*
   * Calculates the structure factor, which is computed as a Fourier Transform.
   * The structure factor is more useful in FCC and BCC lattices
   */

  double s = pow((N / _rho), (1.0 / 3.0));
  double fkx1 = 2.0 * PI / (s / 2.0 * Nx);
  double fky1 = 2.0 * PI / (s / 2.0 * Ny);
  double fkz1 = 2.0 * PI / (s / 2.0 * Nz);
  double sfcosx = 0, sfcosy = 0, sfcosz = 0;
  double sfsinx = 0, sfsiny = 0, sfsinz = 0;

  /* Try to calculate the structure factor at once for all axis
     if the particles per axis are equal. Else simply do them individually */

  if (Nx == Ny && Nx == Nz) {
    for (size_t i = 0; i < rx.size(); ++i) {
      sfcosx += cos(fkx1 * rx[i]);
      sfsinx += sin(fkx1 * rx[i]);
      sfcosy += cos(fky1 * ry[i]);
      sfsiny += sin(fky1 * ry[i]);
      sfcosz += cos(fkz1 * rz[i]);
      sfsinz += sin(fkz1 * rz[i]);
    }
  } else {
    for (const auto &i : rx) {
      sfcosx += cos(fkx1 * i);
      sfsinx += sin(fkx1 * i);
    }
    for (const auto &i : ry) {
      sfcosy += cos(fky1 * i);
      sfsiny += sin(fky1 * i);
    }
    for (const auto &i : ry) {
      sfcosz += cos(fkz1 * i);
      sfsinz += sin(fkz1 * i);
    }
  }

  double kx = (pow(sfcosx, 2) + pow(sfsinx, 2)) / N;
  double ky = (pow(sfcosy, 2) + pow(sfsiny, 2)) / N;
  double kz = (pow(sfcosz, 2) + pow(sfsinz, 2)) / N;

  sfx.push_back(kx);
  sfy.push_back(ky);
  sfz.push_back(kz);
}

// MD Simulation
void MD::Simulation(double DENSITY, double TEMPERATURE, double POWER,
                    double A_CST, std::string pp_type) {
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
   * @param pp_type: The type of the pair potential the simulation is
   *                 modelling. Options are "BIP", "GCM", "EXP"
   */
  // Initialise scaling variables
  std::cout << "***************************\n"
               "** MD Simulation started **\n"
               "***************************\n"
            << std::endl;
  std::cout << "Fluid parameters: rho: " << DENSITY << " T: " << TEMPERATURE
            << " n: " << POWER << " a: " << A_CST << std::endl;

  _rho = DENSITY;
  _T0 = TEMPERATURE;
  dt = 0.005 / sqrt(_T0);  // Time-step, defined here and reused in the Verlet
  // Box length scaling
  scale = pow((N / _rho), (1.0 / 3.0)) / pow(N, (1.0 / 3.0));
  L = pow((N / _rho), 1.0 / 3.0);
  Vol = N / _rho;

  // cut_off redefinition
  // NOTE: Large cut offs increase the runtime exponentially
  cut_off = 3.0;  // TODO: return as argument from BIP or add calibration method
  rg = cut_off;
  dr = rg / nhist;

  // Gets the pair potential for the simulation based on a map of the
  // initials of the pair potential and the pair potential itself.
  pair_potential_type pair_potential_force = get_force_func(pp_type);

  // Generating the filenames for the output
  // Start a new stream only if the fluid is not being compressed
  if (c_counter == 0) {
    std::string data =
        logger.file_naming("/Data", DENSITY, TEMPERATURE, POWER, A_CST);
    std::string pos = logger.file_naming("/Positions_Velocities", DENSITY,
                                         TEMPERATURE, POWER, A_CST);
    std::string rdf =
        logger.file_naming("/RDF", DENSITY, TEMPERATURE, POWER, A_CST);

    logger.open_files(data, rdf, pos);
    logger.time_stamp(
        logger.DATA,
        "# step \t rho \t T \t U \t K \t Pc \t Pk \t MSD \t VAF \t SFx "
        "\t SFy \t SFz");
  }

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);

  for (_STEP_INDEX = 0; _STEP_INDEX < STEPS; ++_STEP_INDEX) {
    // Forces loop
    // Resetting forces
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    U = 0;  // reseting <Potential> U to 0
    PC = 0;

    size_t i, j;
    for (i = 0; i < N - 1; ++i) {
      for (j = i + 1; j < N; ++j) {
        x = rx[i] - rx[j];  // Separation distance
        y = ry[i] - ry[j];  // between particles i and j
        z = rz[i] - rz[j];  // in Cartesian

        // Get the shortest image of the two particles
        // if the particles are near the periodic boundary,
        // this ismage is their reflection.
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
          // Allows the user to choose different pair potentials
          auto [ff, temp_u] = pair_potential_force(r, POWER, A_CST);

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
    U /= N;  // Ensuring that variables are not optimised away
    u_en.push_back(U);

    // Average Configurational Pressure Pc
    PC /= (3 * Vol);  // Ensuring that variables are not optimised away
    pc.push_back(PC);

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

    // Calculate the structure factor k-vectors
    structure_factor(rx, ry, rz);

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
        *pos_x,
        logger.file_naming("/x_data", DENSITY, TEMPERATURE, POWER, A_CST), "\t",
        true);
    f.Write2File<double>(
        *pos_y,
        logger.file_naming("/y_data", DENSITY, TEMPERATURE, POWER, A_CST), "\t",
        true);
    f.Write2File<double>(
        *pos_z,
        logger.file_naming("/z_data", DENSITY, TEMPERATURE, POWER, A_CST), "\t",
        true);
  }

  logger.write_data_file(density, temperature, u_en, k_en, pc, pk, msd, Cr, sfx,
                         sfy, sfz);
  // Saving Last Position
  // todo: saves the same shit
  logger.time_stamp(logger.POS,
                    "# X\tY\tZ\tVx\tVy\tVz\tFx\tFy\tFz");  // fstream
  for (size_t el = 0; el < rx.size(); ++el) {
    logger.POS << rx[el] << '\t' << ry[el] << '\t' << rz[el] << '\t' << vx[el]
               << '\t' << vy[el] << '\t' << vz[el] << '\t' << fx[el] << '\t'
               << fy[el] << '\t' << fz[el] << std::endl;
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

  // Close file streams, makes Simulation reusable in loops
  reset_values();
}

// todo: move to other class devoted to phase transitions and multiphase
void MD::get_phases(double DENSITY, double FINAL_DENSITY, double DENSITY_INC,
                    double TEMPERATURE, double POWER, double A_CST,
                    std::string pp_type) {
  /*
   * Compress the fluid to get the phase boundary for a specific temperature.
   *
   * Performs repeated compresss of the fluid by periodically
   * incrementing the density of the fluid.
   * As a consequence the box length, the scaling factor and the
   * position vectors are also scaled in order to conserve the number
   * of particles in the box.
   *
   */
  // todo: many features do not work at this moment like particle tracking
  // todo: the POS << time_stamp stream will be a mess, same with RDF
  double current_rho = DENSITY;
  double old_box_length = 0;
  // size_t compression_num = ceil((FINAL_DENSITY - DENSITY) / DENSITY_INC);

  do {
    Simulation(current_rho, TEMPERATURE, POWER, A_CST, pp_type);
    // Holds the box length of the previous simulation just run
    old_box_length = L;
    // Density incrementation
    current_rho += DENSITY_INC;
    // Simulation updates old_box_length
    // the updated current_rho can generate the new box length
    // This value gets recalculated in the next Simulation
    L = pow((N / current_rho), 1.0 / 3.0);
    double box_length_ratio = L / old_box_length;

    // Rescalling the positional vectors
    for (size_t i = 0; i < N; ++i) {
      rx[i] *= box_length_ratio;
      ry[i] *= box_length_ratio;
      rz[i] *= box_length_ratio;
    }
    ++c_counter;
  } while (abs(current_rho - FINAL_DENSITY) > 0.00001);
  reset_values();
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

  /*
    todo: add try-except block for stream closing,
    bc the stream might be closed and the user might then call reset_values
    which will throw an exception
  */

  // Do not close streams and do not clear position and velocity vectors
  // in the case where the fluid is being compressed
  if (!compress) {
    // Close streams
    logger.RDF.close();
    logger.DATA.close();
    logger.POS.close();
    // Clear values, size, but reserve capacity
    rx.clear();
    ry.clear();
    rz.clear();
    vx.clear();
    vy.clear();
    vz.clear();
  }
  // Reset the MSD initial vectors
  rrx.clear();
  rry.clear();
  rrz.clear();
  // Clear monitored quantities
  density.clear();
  temperature.clear();
  u_en.clear();
  k_en.clear();
  pc.clear();
  pk.clear();
  msd.clear();
  Cr.clear();
  gr.resize(nhist + 1, 0);  // gr with Index igr
}
