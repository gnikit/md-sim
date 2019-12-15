#include "MD.h"
// TODO: create internal structures for quantities such as velocities, rs, fs,
// todo: add logger https://github.com/gabime/spdlog

// FileIO has to be loaded after the math libraries
#include "FileIO.h"  // FileIO class

#pragma warning(disable : 4996)  //_CRT_SECURE_NO_WARNINGS

MD::MD(size_t step_number, std::vector<size_t> particles, std::string lattice) {
  // Assign number of iterations of the MD algorithm
  __steps = step_number;
  std::cout << "Number of steps: " << __steps << std::endl;

  // Assign the type of lattice
  __lattice = lattice;
  std::cout << "Lattice type: " << __lattice << std::endl;

  try {
    if (particles.empty()) {
      throw "The supplied particles vector is empty";
    } else if (particles.size() < 3) {
      throw "The supplied particles vector is of incorrect size";
    } else if (std::find(particles.begin(), particles.end(), 0) !=
               particles.end()) {
      throw
                  "The supplied particles vector contains a 0\n"
                  "particles cannot be 0 in x, y or z";
    }

  } catch (const char *msg) {
    std::cerr << "Error: " << msg << std::endl;
    exit(1);
  }

  // Calculate the total number of particles N based on the lattice
  __Nx = particles[0];
  __Ny = particles[1];
  __Nz = particles[2];
  if (lattice == "FCC") {
    __N = __Nx * __Ny * __Nz * 4;
  } else if (lattice == "BCC") {
    __N = __Nx * __Ny * __Nz * 2;
  } else {
    __N = __Nx * __Ny * __Nz;
  }
  std::cout << "Number of particles: " << __N << std::endl;

  __dir = ".";

  // If compress is true, then STEPS = steps_per_compression
  __compress = false;

  // Save all the positions for the fluid
  __visualise = false;

  // Accuracy of RDF
  __nhist = 500;

  __rdf_wait = 0;

  // For efficiency, memory in the containers is reserved before use
  /* Positions */
  rx.reserve(__N);
  ry.reserve(__N);
  rz.reserve(__N);
  /* Velocities */
  vx.reserve(__N);
  vy.reserve(__N);
  vz.reserve(__N);
  /* RDF */
  gr.resize(__nhist + 1, 0);  // gr with Index igr
  /* Forces/Acceleration */
  fx.resize(__N, 0);
  fy.resize(__N, 0);
  fz.resize(__N, 0);
  /* Structure factor k-arrays */
  sfx.reserve(__N);
  sfy.reserve(__N);
  sfz.reserve(__N);
  /* Observed Quantities */
  Cr.reserve(__steps);    // Velocity Autocorrelation Function
  msd.reserve(__steps);   // Mean Square Displacement
  u_en.reserve(__steps);  // Average Potential Energy
  k_en.reserve(__steps);  // Average Kinetic Energy
  pc.reserve(__steps);    // Configurational Pressure
  pk.reserve(__steps);    // Kinetic Pressure
  temperature.reserve(__steps);

  /* Visualisation vectors on the heap*/
  pos_x = new std::vector<std::vector<double>>(__steps);
  pos_y = new std::vector<std::vector<double>>(__steps);
  pos_z = new std::vector<std::vector<double>>(__steps);

  PI = acos(-1.0);
  fixed_seed = false;
}

// Delegating constructors with reduced number of arguments
// https://en.wikipedia.org/wiki/C++11#Object_construction_improvement
// Convinient constructor to use for simple cases

MD::~MD() {
  // Destroy the vectors allocated on the heap
  delete pos_x;
  delete pos_y;
  delete pos_z;
}

void MD::load_options(std::string out_directory = ".",
                      bool track_particles = false, size_t rdf_bins = 500,
                      size_t collect_rdf_after = 500) {
  // Test whether the input directory exists
  if (!out_directory.empty()) {
    try {
      __dir = out_directory;
      if (!fs::exists(__dir)) {
        throw
          "input out_directory in MD constructor does not exist.\n"
          "Use a valid directory for output files to be saved";
      }

    } catch (const char *msg) {
      std::cerr << "Error: " << msg << std::endl;
      exit(1);
    }
  }

  std::cout << "Output directory set to: " << __dir << std::endl;

  // Save all the positions for the fluid
  __visualise = track_particles;
  std::cout << "Particle visualisation set to: " << __visualise << std::endl;

  // Accuracy of RDF
  __nhist = rdf_bins;
  std::cout << "RDF accuracy set to: " << __nhist << " bins" << std::endl;

  // Ensuring the number of steps is greater than the rdf equilibration period
  try {
    /* The number of iterations the data collection of RDF is postponed
       in order to allow the fluid to lose its internal cubic lattice */
    __rdf_wait = collect_rdf_after;

    /* Substraction of size_ts if negative results into logic errors
       hence the use of an int temp; */
    int temp = __steps - collect_rdf_after;
    if (temp < 0) {
      throw "collect_rdf_after is greater than the step_number";
    }
  } catch (const char *msg) {
    std::cerr << "Warning: " << msg << std::endl;
    std::cerr << "         rdf_wait is set to 0" << std::endl;
    __rdf_wait = 0;
  }
  std::cout << "RDF equilibration period set to: " << __rdf_wait << std::endl;
}

// Methods for MD Analysis
double MD::initialise(std::vector<double> &x, std::vector<double> &y,
                      std::vector<double> &z, std::vector<double> &vx,
                      std::vector<double> &vy, std::vector<double> &vz,
                      double TEMPERATURE) {
  // Initialise position matrix and velocity matrix from Cubic Centred Lattice
  if (!__compress || (__compress && c_counter == 0)) {
    if (__lattice == "FCC") {
      // Coordinates for the FCC lattice
      double x_c[4] = {0.25, 0.75, 0.75, 0.25};
      double y_c[4] = {0.25, 0.75, 0.25, 0.75};
      double z_c[4] = {0.25, 0.25, 0.75, 0.75};

      // Loop over the the corner coordinates of the FCC and then x, y, z
      for (size_t c = 0; c < 4; ++c) {
        for (size_t i = 0; i < __Nx; ++i) {
          for (size_t j = 0; j < __Ny; ++j) {
            for (size_t k = 0; k < __Nz; ++k) {
              x.push_back((i + x_c[c]) * (__L / __Nx));
              y.push_back((j + y_c[c]) * (__L / __Ny));
              z.push_back((k + z_c[c]) * (__L / __Nz));
            }
          }
        }
      }
    }

    // BCC lattice
    else if (__lattice == "BCC") {
      double x_c[2] = {0.25, 0.75};
      double y_c[2] = {0.25, 0.75};
      double z_c[2] = {0.25, 0.75};

      for (size_t c = 0; c < 2; ++c) {
        for (size_t i = 0; i < __Nx; i++) {
          for (size_t j = 0; j < __Ny; j++) {
            for (size_t k = 0; k < __Nz; k++) {
              x.push_back((i + x_c[c]) * (__L / __Nx));
              y.push_back((j + y_c[c]) * (__L / __Ny));
              z.push_back((k + z_c[c]) * (__L / __Nz));
            }
          }
        }
      }
    }

    else if (__lattice == "RANDOM") {
      x.resize(__N);
      y.resize(__N);
      z.resize(__N);
      // todo: add schema option to set variance
      mb_distribution(x, y, z, TEMPERATURE);
    }

    // Simple Cubic lattice
    else {
      for (size_t i = 0; i < __Nx; ++i) {
        for (size_t j = 0; j < __Ny; ++j) {
          for (size_t k = 0; k < __Nz; ++k) {
            x.push_back((i + 0.5) * (__L / __Nx));
            y.push_back((j + 0.5) * (__L / __Ny));
            z.push_back((k + 0.5) * (__L / __Nz));
          }
        }
      }
    }

    // Generates Maxwell-Boltzmann distribution
    mb_distribution(vx, vy, vz, TEMPERATURE);
  }

  // Calculate the average velocities
  double mean_vx = std::accumulate(vx.begin(), vx.end(), 0.0) / __N;
  double mean_vy = std::accumulate(vy.begin(), vy.end(), 0.0) / __N;
  double mean_vz = std::accumulate(vz.begin(), vz.end(), 0.0) / __N;
  // Conserve the momentum of the fluid by subsracting the average velocities
  // using a lambda expression
  std::for_each(vx.begin(), vx.end(), [mean_vx](double &d) { d -= mean_vx; });
  std::for_each(vy.begin(), vy.end(), [mean_vy](double &d) { d -= mean_vy; });
  std::for_each(vz.begin(), vz.end(), [mean_vz](double &d) { d -= mean_vz; });

  size_t i;
  // Temperature calculation, statistically
  double KE = 0;
  for (i = 0; i < __N; ++i) {
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  __T = KE / (1.5 * __N);
  scale_v = sqrt(TEMPERATURE / __T);  // scaling factor

  // Velocity scaling
  for (i = 0; i < __N; ++i) {
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

  return KE;
}

void MD::mb_distribution(std::vector<double> &vx, std::vector<double> &vy,
                         std::vector<double> &vz, double TEMPERATURE) {
  double kb = 1.0;
  double m = 1.0;

  double var = sqrt(TEMPERATURE * kb / m);
  double mean = 0;

  // Use current time as seed for random generator
  std::srand(std::time(nullptr));
  int random_variable = std::rand();
  if (fixed_seed) random_variable = 666;  // Fixing it for testing

  std::default_random_engine generator;
  generator.seed(random_variable);

  std::normal_distribution<double> g_x(mean, var);
  std::normal_distribution<double> g_y(mean, var);
  std::normal_distribution<double> g_z(mean, var);

  for (size_t i = 0; i < __N; ++i) {
    vx.push_back(g_x(generator));
    vy.push_back(g_y(generator));
    vz.push_back(g_z(generator));
  }
}

void MD::set_random_position_variance(double var) { __var = var * var; }

double MD::verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry,
                            std::vector<double> &rz, std::vector<double> &vx,
                            std::vector<double> &vy, std::vector<double> &vz,
                            bool sample_msd = true) {
  size_t i;
  double KE = 0;

  for (i = 0; i < __N; ++i) {
    // Step velocities forward in time
    vx[i] = vx[i] * scale_v + fx[i] * __dt;
    vy[i] = vy[i] * scale_v + fy[i] * __dt;
    vz[i] = vz[i] * scale_v + fz[i] * __dt;

    // Step positions forward in time
    rx[i] = rx[i] + vx[i] * __dt;
    ry[i] = ry[i] + vy[i] * __dt;
    rz[i] = rz[i] + vz[i] * __dt;

    if (sample_msd) {
      // MSD stepping
      rrx[i] = rrx[i] + vx[i] * __dt;
      rry[i] = rry[i] + vy[i] * __dt;
      rrz[i] = rrz[i] + vz[i] * __dt;
    }

    // Kinetic Energy Calculation
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    // Apply periodic boundary conditions to ensure particles remain
    // inside the box
    // todo: make boundary conditions routines
    if (rx[i] > __L) rx[i] = rx[i] - __L;
    if (rx[i] < 0.0) rx[i] = rx[i] + __L;
    if (ry[i] > __L) ry[i] = ry[i] - __L;
    if (ry[i] < 0.0) ry[i] = ry[i] + __L;
    if (rz[i] > __L) rz[i] = rz[i] - __L;
    if (rz[i] < 0.0) rz[i] = rz[i] + __L;
  }

  return KE;
}

void MD::velocity_autocorrelation_function(std::vector<double> &Cvx,
                                           std::vector<double> &Cvy,
                                           std::vector<double> &Cvz,
                                           std::vector<double> &vx,
                                           std::vector<double> &vy,
                                           std::vector<double> &vz) {
  double cr_temp = 0;  // resets the sum every time step
  double m = 1.0;      // particle mass
  size_t i;
  /* The peak of the VAF is located at 3kb*T/m */
  double norm = 3 * __T0 / m;
  for (i = 0; i < __N; i++) {
    cr_temp += (Cvx[i] * vx[i] + Cvy[i] * vy[i] + Cvz[i] * vz[i]);
  }
  Cr.push_back((cr_temp / __N) / norm);
}

void MD::radial_distribution_function(double &rho, double &cut_off,
                                      size_t &bins, size_t &particles) {
  double R = 0;
  double norm = 1;
  // Exclude the self particle interaction from the density
  double cor_rho = rho * (particles - 1) / particles;
  double dr = cut_off / bins;

  logger.RDF << "# particles (N): " << particles << " steps: " << __steps
             << " rho: " << rho << " bins: " << bins
             << " cut_off (rg): " << cut_off << " dr: " << dr << std::endl;
  logger.RDF << "# Radius (r)" << '\t' << "Normalised" << '\t' << "Unormalised"
             << std::endl;

  for (size_t i = 1; i < bins; ++i) {
    R = cut_off * i / bins;
    // Volume between 2 spheres, accounting for double counting
    // hence the 2/3*pi*((R+dr)**3 - R**3)
    // Accounting for the rdf_wait time steps
    norm = cor_rho * (2.0 / 3.0 * PI * particles * (__steps - __rdf_wait) *
                      (pow((R + (dr / 2.0)), 3) - pow((R - (dr / 2.0)), 3)));

    logger.RDF << R << '\t' << gr[i] / norm << '\t' << gr[i] << std::endl;
  }
}

void MD::mean_square_displacement(std::vector<double> &MSDx,
                                  std::vector<double> &MSDy,
                                  std::vector<double> &MSDz,
                                  std::vector<double> &rrx,
                                  std::vector<double> &rry,
                                  std::vector<double> &rrz) {
  double msd_temp = 0;

  for (size_t i = 0; i < __N; ++i) {
    msd_temp += (pow((rrx[i] - MSDx[i]), 2) + pow((rry[i] - MSDy[i]), 2) +
                 pow((rrz[i] - MSDz[i]), 2));
  }
  msd.push_back(msd_temp / __N);
}

void MD::structure_factor(std::vector<double> &rx, std::vector<double> &ry,
                          std::vector<double> &rz) {
  double s = pow((__N / __rho), (1.0 / 3.0));
  double fkx1 = 2.0 * PI / (s / (2.0 * __Nx));
  double fky1 = 2.0 * PI / (s / (2.0 * __Ny));
  double fkz1 = 2.0 * PI / (s / (2.0 * __Nz));
  double sfcosx = 0, sfcosy = 0, sfcosz = 0;
  double sfsinx = 0, sfsiny = 0, sfsinz = 0;

  /* Try to calculate the structure factor at once for all axis
     if the particles per axis are equal. Else simply do them individually */

  if (__Nx == __Ny && __Nx == __Nz) {
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

  double kx = sqrt(pow(sfcosx / __N, 2) + pow(sfsinx / __N, 2));
  double ky = sqrt(pow(sfcosy / __N, 2) + pow(sfsiny / __N, 2));
  double kz = sqrt(pow(sfcosz / __N, 2) + pow(sfsinz / __N, 2));

  sfx.push_back(kx);
  sfy.push_back(ky);
  sfz.push_back(kz);
}

void MD::simulation(std::string simulation_name, double DENSITY,
                    double TEMPERATURE, double POWER = NAN, double A_CST = NAN,
                    std::string pp_type = "LennardJones") {
  // Initialise the variables with the input parameters
  // Name the simulation. This will be used as a prefix in the files
  __simulation_name = simulation_name;
  __rho = DENSITY;
  __T0 = TEMPERATURE;
  __power = POWER;
  __a_cst = A_CST;
  __pp_type = pp_type;

  std::cout << "***************************\n"
               "** MD simulation started **\n"
               "***************************\n"
            << std::endl;

  // Sets the unneeded variables (A and/or n) to NAN depending on the pp-type
  std::cout << set_simulation_params(DENSITY, TEMPERATURE, POWER, A_CST,
                                     pp_type)
            << std::endl;

  // Initialise scaling variables
  __dt = 0.005 / sqrt(__T0);  // dt defined here and reused in the Verlet
  // Box length scaling
  __L = pow((__N / __rho), 1.0 / 3.0);
  double Vol = __N / __rho;

  // cut_off definition
  // Hard coded into 1/3 of the box length
  // NOTE: Large cut offs increase the runtime exponentially
  __cut_off = __L / 3.0;  // TODO: return as arg from pp or add calibration func
  // if cut-off is too large rescale it
  if (__cut_off > __L / 2.0) {
    std::cerr << "Warning: cutoff was too large!\n"
                 "Setting cut-off to half the length box\n"
                 "cut-off: "
              << __L / 2.0 << std::endl;
    __cut_off = __L / 2.0;
  }

  /* Gets the pair potential for the simulation based on a map of the
     initials of the pair potential and the pair potential itself. */
  pair_potential_type pair_potential_force = get_force_func(pp_type);

  // Generating the filenames for the output
  // Start a new stream only if the fluid is not being compressed
  if (c_counter == 0) {
    std::string data =
        __dir + logger.file_naming("/" + __simulation_name + "Data", __steps,
                                   __N, __rho, __T0, __power, __a_cst);
    std::string pos =
        __dir +
        logger.file_naming("/" + __simulation_name + "Positions_Velocities",
                           __steps, __N, __rho, __T0, __power, __a_cst);
    std::string rdf =
        __dir + logger.file_naming("/" + __simulation_name + "RDF", __steps,
                                   __N, __rho, __T0, __power, __a_cst);

    logger.open_files(data, rdf, pos); //todo: make more general
    logger.time_stamp(logger.DATA,
                      "# step \t rho \t T \t U \t K \t Pc \t Pk \t MSD \t VAF "
                      "\t SFx \t SFy \t SFz");
  }

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  // Initialise the simulation, lattice params and much more
  __KE = initialise(rx, ry, rz, vx, vy, vz, __T0);

  for (__step_idx = 0; __step_idx < __steps; ++__step_idx) {
    // Forces loop
    // Resetting forces
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    // Reseting <Potential> U to 0
    double U = 0;   // Potential Energy
    double PC = 0;  // Configurational Pressure

    size_t i, j;
    for (i = 0; i < __N - 1; ++i) {
      for (j = i + 1; j < __N; ++j) {
        // distance between particle i and j
        double x = rx[i] - rx[j];  // Separation distance
        double y = ry[i] - ry[j];  // between particles i and j
        double z = rz[i] - rz[j];  // in Cartesian

        // Get the shortest image of the two particles
        // if the particles are near the periodic boundary,
        // this image is their reflection.
        if (x > (0.5 * __L)) x = x - __L;
        if (x < (-0.5 * __L)) x = x + __L;
        if (y > (0.5 * __L)) y = y - __L;
        if (y < (-0.5 * __L)) y = y + __L;
        if (z > (0.5 * __L)) z = z - __L;
        if (z < (-0.5 * __L)) z = z + __L;

        // Pair potential radius
        double r = sqrt((x * x) + (y * y) + (z * z));

        // Force loop
        if (r < __cut_off) {
          // Allows the user to choose different pair potentials
          auto [ff, temp_u] = pair_potential_force(r, __power, __a_cst);

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
          if (__step_idx > __rdf_wait) {
            igr = round(__nhist * r / __cut_off);
            gr[igr] += 1;
          }
        }
      }
    }

    // Average Potential Energy per particle
    u_en.push_back(U / __N);

    // Average Configurational Pressure Pc
    pc.push_back(PC / (3 * Vol));

    // Isothermal Calibration
    scale_v = sqrt(__T0 / __T);  // using T & KE from prev timestep

    __KE = verlet_algorithm(rx, ry, rz, vx, vy, vz, true);
    mean_square_displacement(MSDx, MSDy, MSDz, rrx, rry, rrz);
    velocity_autocorrelation_function(Cvx, Cvy, Cvz, vx, vy, vz);

    // Average Temperature
    __T = __KE / (1.5 * __N);
    temperature.push_back(__T);

    // Kinetic Pressure
    pk.push_back(__rho * __T);

    // Average Kinetic Energy
    k_en.push_back(__KE / __N);

    // Density
    density.push_back(__rho);

    // Calculate the structure factor k-vectors
    structure_factor(rx, ry, rz);

    // Save positions for visualisation with Python
    if (__visualise) {
      // Reserve memory for the position vectors
      (*pos_x)[__step_idx].reserve(__N);
      (*pos_y)[__step_idx].reserve(__N);
      (*pos_z)[__step_idx].reserve(__N);

      // Populate the vectors with the current positions
      (*pos_x)[__step_idx] = rx;
      (*pos_y)[__step_idx] = ry;
      (*pos_z)[__step_idx] = rz;
    }
  }
  // simulation Ends HERE

  if (__visualise) {
    // Save particle positions to files
    FileIO f;
    // Write the arrays as jagged,(hence transposed), this creates rows=STEPS
    // and columns=PARTICLES
    f.Write2File<double>(
        *pos_x,
        logger.file_naming(__dir + "/" + __simulation_name + "x_data", __steps,
                           __N, __rho, __T0, __power, __a_cst),
        "\t", true);
    f.Write2File<double>(
        *pos_y,
        logger.file_naming(__dir + "/" + __simulation_name + "y_data", __steps,
                           __N, __rho, __T0, __power, __a_cst),
        "\t", true);
    f.Write2File<double>(
        *pos_z,
        logger.file_naming(__dir + "/" + __simulation_name + "z_data", __steps,
                           __N, __rho, __T0, __power, __a_cst),
        "\t", true);
  }

  logger.write_data_file(__steps, density, temperature, u_en, k_en, pc, pk, msd,
                         Cr, sfx, sfy, sfz);
  // Saving Last Position
  // todo: if we are compressing, save the last position of the compression step
  logger.time_stamp(logger.POS, "# X\tY\tZ\tVx\tVy\tVz\tFx\tFy\tFz");

  for (size_t el = 0; el < rx.size(); ++el) {
    logger.POS << rx[el] << '\t' << ry[el] << '\t' << rz[el] << '\t' << vx[el]
               << '\t' << vy[el] << '\t' << vz[el] << '\t' << fx[el] << '\t'
               << fy[el] << '\t' << fz[el] << std::endl;
  }

  radial_distribution_function(__rho, __cut_off, __nhist, __N);

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
               "** MD simulation terminated **\n"
               "******************************\n"
            << std::endl;

  // Close file streams, makes simulation reusable in loops
  reset_values();
}

void MD::reset_values(bool force_reset) {
  /*
    bc the stream might be closed and the user might then call reset_values
    which will throw an exception
  */

  // Do not close streams and do not clear position and velocity vectors
  // in the case where the fluid is being compressed
  if (!__compress || force_reset) {
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
    c_counter = 0;
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
  gr.resize(__nhist + 1, 0);  // gr with Index igr
}

std::string MD::get_dir() { return __dir; }

std::string MD::get_simulation_name() { return __simulation_name; }

std::string MD::set_simulation_params(double &rho, double &T, double &power,
                                      double &a, std::string &pp_type) {
  std::string params =
      "Fluid parameters: rho: " + stat_file::convert_to_string(rho, 4) +
      " T: " + stat_file::convert_to_string(T, 4);

  params = "Lattice: " + get_lattice_structure() + "\n" + params;

  if (pp_type == "GaussianCoreModel") {
    params = "Potential: GaussianCoreModel, " + params;
    __power = NAN;  // Set the variable to NAN to be ignore by the logger
    __a_cst = NAN;  // Set the variable to NAN to be ignore by the logger
  }

  else if (pp_type == "LennardJones") {
    params = "Potential: LennardJones, " + params;
    __power = NAN;  // Set the variable to NAN to be ignore by the logger
    __a_cst = NAN;  // Set the variable to NAN to be ignore by the logger
  }

  else if (pp_type == "Exponential") {
    params = "Potential: Exponential, " + params;
    params += " m: " + stat_file::convert_to_string(power, 4);
    params += " C: " + stat_file::convert_to_string(a, 4);
  }

  else if (pp_type == "BoundedInversePower") {
    params = "Potential: BoundedInversePower, " + params;
    params += " n: " + stat_file::convert_to_string(power, 4);
    params += " A: " + stat_file::convert_to_string(a, 4);
  }

  else {
    std::cerr << "Warning pp_type unknown\n"
              << "Defaulting to BIP potential" << std::endl;
    params = "Potential: BIP, " + params;
    params += " n: " + stat_file::convert_to_string(power, 4);
    params += " A: " + stat_file::convert_to_string(a, 4);
  }

  return params;
}

std::string MD::get_lattice_structure() { return __lattice; }

bool MD::get_visualisation_flag() { return __visualise; }

void MD::set_visualisation_flag(bool is_visualising) {
  __visualise = is_visualising;
}

size_t MD::get_particle_number() { return __N; }

size_t MD::get_rdf_accuracy() { return __nhist; }

void MD::enable_testing(bool is_testing) { fixed_seed = is_testing; }

size_t MD::get_rdf_collect_after() { return __rdf_wait; }

void MD::set_rdf_collect_after(size_t rdf_collect_after) {
  __rdf_wait = rdf_collect_after;
}