#include "MD.h"
// TODO: create internal structures for quantities such as velocities, rs, fs,
// todo: add logger https://github.com/gabime/spdlog
// TODO: scale the box by Lx, Ly, Lz in a tensor form

MD::MD() {}

MD::MD(options_type &input_options) {
  /* Test whether the input directory exists */
  if (!input_options.io_options.dir.empty()) {
    try {
      options.io_options.dir = input_options.io_options.dir;
      if (!fs::exists(options.io_options.dir)) {
        throw
          "input out_directory in MD constructor does not exist.\n"
          "Use a valid directory for output files to be saved";
      }

    } catch (const char *msg) {
      std::cerr << "Error: " << msg << std::endl;
      exit(1);
    }
  }
  std::cout << "Output directory set to: " << options.io_options.dir
            << std::endl;

  /* Pass type of simulation */
  options.simulation_type = input_options.simulation_type;
  std::cout << "Simulation type: " << options.simulation_type << std::endl;

  /* Pass simulation name if any */
  options.io_options.simulation_name = input_options.io_options.simulation_name;
  std::cout << "Simulation name: " << options.io_options.simulation_name
            << std::endl;

  /* Save all the positions for the fluid */
  options.io_options.visualise = input_options.io_options.visualise;
  std::cout << "Particle visualisation: " << options.io_options.visualise
            << std::endl;

  /* Pass stepping algorithm */
  // todo: test string against availbale options
  options.stepping_alg = input_options.stepping_alg;
  std::cout << "Iterative algorithm: " << options.stepping_alg << std::endl;

  /* Pass number of iterations */
  options.steps = input_options.steps;
  std::cout << "Number of steps: " << options.steps << std::endl;

  /* Pass particles and lattice */
  options.lattice = input_options.lattice;
  std::cout << "Initial lattice: " << options.lattice << std::endl;

  /* Pass particles */
  try {
    if (input_options.particles.empty()) {
      throw "The supplied particles vector is empty";
    } else if (input_options.particles.size() < 3) {
      throw "The supplied particles vector is of incorrect size";
    } else if (std::find(input_options.particles.begin(),
                         input_options.particles.end(),
                         0) != input_options.particles.end()) {
      throw
            "The supplied particles vector contains a 0\n"
            "particles cannot be 0 in x, y or z";
    }

  } catch (const char *msg) {
    std::cerr << "Error: " << msg << std::endl;
    exit(1);
  }
  options.particles = input_options.particles;

  /* Calculate the total number of particles N based on the lattice */
  options.Nx = input_options.particles[0];
  options.Ny = input_options.particles[1];
  options.Nz = input_options.particles[2];
  if (input_options.lattice == "FCC") {
    options.N = options.Nx * options.Ny * options.Nz * 4;
  } else if (input_options.lattice == "BCC") {
    options.N = options.Nx * options.Ny * options.Nz * 2;
  } else {
    options.N = options.Nx * options.Ny * options.Nz;
  }
  std::cout << "Number of particles: " << options.N << std::endl;

  /* Pass physical parameters */
  /* Pass the pair potential */
  options.potential_type = input_options.potential_type;
  std::cout << "Pair potential: " << options.potential_type << std::endl;

  if (input_options.density > 0)
    options.density = input_options.density;
  else {
    std::cerr << "Error: Negative density supplied" << std::endl;
    exit(1);
  }

  if (input_options.target_temperature > 0)
    options.target_temperature = input_options.target_temperature;
  else {
    std::cerr << "Error: Negative temperature supplied" << std::endl;
    exit(1);
  }
  options.power = input_options.power;

  options.a_cst = input_options.a_cst;

  /* Initialise scaling variables */
  options.dt =
      0.005 / sqrt(options.target_temperature); /* todo: add to schema */
  /* Box length scaling */
  options.L = pow((options.N / options.density), 1.0 / 3.0);
  options.Lx = options.Ly = options.Lz = options.L; /* todo: questionable! */
  options.volume = options.N / options.density;

  /* cut_off definition */
  if (input_options.cut_off > 0) {
    options.cut_off = input_options.cut_off;
    /* if cut-off is too large rescale it */
    if (options.cut_off > options.L / 2.0) {
      std::cerr << "Warning: cutoff was too large!\n"
                   "Setting cut-off to half the length box\n"
                   "cut-off: "
                << options.L / 2.0 << std::endl;
      options.cut_off = options.L / 2.0;
    }
  } else {
    /* Hard coded into 1/3 of the box length */
    /* NOTE: Large cut offs increase the runtime exponentially */
    options.cut_off = options.L / 3.0;
  }

  /* Set boundary conditions //todo */

  /* Accuracy of RDF */
  options.rdf_options.rdf_bins = input_options.rdf_options.rdf_bins;
  std::cout << "RDF bins: " << options.rdf_options.rdf_bins << std::endl;

  /* Ensuring the number of steps is greater than the rdf equilibration period
   */
  try {
    /* The number of iterations the data collection of RDF is postponed
       in order to allow the fluid to lose its internal cubic lattice */
    options.rdf_options.rdf_wait = input_options.rdf_options.rdf_wait;

    /* Substraction of size_ts if negative results into logic errors
       hence the use of an int temp; */
    int temp = options.steps - options.rdf_options.rdf_wait;
    if (temp < 0) {
      throw "collect_rdf_after is greater than the step_number";
    }
  } catch (const char *msg) {
    std::cerr << "Warning: " << msg << std::endl;
    std::cerr << "         rdf_wait is set to 0" << std::endl;
    options.rdf_options.rdf_wait = 0;
  }
  std::cout << "RDF equilibration period set to: "
            << options.rdf_options.rdf_wait << std::endl;
  /* Pass testing options */
  options.test_options.is_testing = input_options.test_options.is_testing;
  std::cout << "Testing: " << options.test_options.is_testing << std::endl;

  /* Visualisation vectors on the heap*/
  pos_x = new std::vector<std::vector<double>>(options.steps);
  pos_y = new std::vector<std::vector<double>>(options.steps);
  pos_z = new std::vector<std::vector<double>>(options.steps);

  PI = acos(-1.0);
}

MD::MD(size_t step_number, std::vector<size_t> particles, std::string lattice) {
  /* Assign number of iterations of the MD algorithm */
  options.steps = step_number;
  std::cout << "Number of steps: " << options.steps << std::endl;

  /* Assign the type of lattice */
  options.lattice = lattice;
  std::cout << "Lattice type: " << options.lattice << std::endl;

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

  /* Calculate the total number of particles N based on the lattice */
  options.Nx = particles[0];
  options.Ny = particles[1];
  options.Nz = particles[2];
  if (lattice == "FCC") {
    options.N = options.Nx * options.Ny * options.Nz * 4;
  } else if (lattice == "BCC") {
    options.N = options.Nx * options.Ny * options.Nz * 2;
  } else {
    options.N = options.Nx * options.Ny * options.Nz;
  }
  std::cout << "Number of particles: " << options.N << std::endl;

  options.io_options.dir = ".";

  /* If compress is true, then STEPS = steps_per_compression */
  options.compression_options.compression = false;

  /* Save all the positions for the fluid */
  options.io_options.visualise = false;

  /* Accuracy of RDF */
  options.rdf_options.rdf_bins = 500;

  options.rdf_options.rdf_wait = 0;

  /* For efficiency, memory in the containers is reserved before use */
  /* Positions */
  r.x.reserve(options.N);
  r.y.reserve(options.N);
  r.z.reserve(options.N);
  /* Velocities */
  v.x.reserve(options.N);
  v.y.reserve(options.N);
  v.z.reserve(options.N);
  /* RDF */
  rdf.resize(options.rdf_options.rdf_bins + 1, 0); /* gr with Index igr */
  /* Forces/Acceleration */
  f.x.resize(options.N, 0);
  f.y.resize(options.N, 0);
  f.z.resize(options.N, 0);
  /* Structure factor k-arrays */
  sf.x.reserve(options.N);
  sf.y.reserve(options.N);
  sf.z.reserve(options.N);
  /* Observed Quantities */
  Cr.reserve(options.steps);   /* Velocity Autocorrelation Function */
  msd.reserve(options.steps);  /* Mean Square Displacement */
  u_en.reserve(options.steps); /* Average Potential Energy */
  k_en.reserve(options.steps); /* Average Kinetic Energy */
  pc.reserve(options.steps);   /* Configurational Pressure */
  pk.reserve(options.steps);   /* Kinetic Pressure */
  temperature.reserve(options.steps);

  /* Visualisation vectors on the heap*/
  pos_x = new std::vector<std::vector<double>>(options.steps);
  pos_y = new std::vector<std::vector<double>>(options.steps);
  pos_z = new std::vector<std::vector<double>>(options.steps);

  PI = acos(-1.0);
  options.test_options.is_testing = false;
}

/* Delegating constructors with reduced number of arguments
 https://en.wikipedia.org/wiki/C++11#Object_construction_improvement
 Convinient constructor to use for simple cases
 */

MD::~MD() {
  /* Destroy the vectors allocated on the heap */
  delete pos_x;
  delete pos_y;
  delete pos_z;
}

void MD::choose_lattice_formation(std::string &lattice, vector_3d &r) {
  if (lattice == "FCC") {
    /* Coordinates for the FCC lattice */
    double x_c[4] = {0.25, 0.75, 0.75, 0.25};
    double y_c[4] = {0.25, 0.75, 0.25, 0.75};
    double z_c[4] = {0.25, 0.25, 0.75, 0.75};

    /* Loop over the the corner coordinates of the FCC and then x, y, z */
    for (size_t c = 0; c < 4; ++c) {
      for (size_t i = 0; i < options.Nx; ++i) {
        for (size_t j = 0; j < options.Ny; ++j) {
          for (size_t k = 0; k < options.Nz; ++k) {
            r.x.push_back((i + x_c[c]) * (options.Lx / options.Nx));
            r.y.push_back((j + y_c[c]) * (options.Ly / options.Ny));
            r.z.push_back((k + z_c[c]) * (options.Lz / options.Nz));
          }
        }
      }
    }
  }

  /* BCC lattice */
  else if (lattice == "BCC") {
    double x_c[2] = {0.25, 0.75};
    double y_c[2] = {0.25, 0.75};
    double z_c[2] = {0.25, 0.75};

    for (size_t c = 0; c < 2; ++c) {
      for (size_t i = 0; i < options.Nx; i++) {
        for (size_t j = 0; j < options.Ny; j++) {
          for (size_t k = 0; k < options.Nz; k++) {
            r.x.push_back((i + x_c[c]) * (options.Lx / options.Nx));
            r.y.push_back((j + y_c[c]) * (options.Ly / options.Ny));
            r.z.push_back((k + z_c[c]) * (options.Lz / options.Nz));
          }
        }
      }
    }
  }

  else if (lattice == "RANDOM") {
    r.x.resize(options.N);
    r.y.resize(options.N);
    r.z.resize(options.N);
    // todo: test
    mb_distribution(r, options.random_lattice_var);
  }

  /* Simple Cubic lattice */
  else {
    for (size_t i = 0; i < options.Nx; ++i) {
      for (size_t j = 0; j < options.Ny; ++j) {
        for (size_t k = 0; k < options.Nz; ++k) {
          r.x.push_back((i + 0.5) * (options.Lx / options.Nx));
          r.y.push_back((j + 0.5) * (options.Ly / options.Ny));
          r.z.push_back((k + 0.5) * (options.Lz / options.Nz));
        }
      }
    }
  }
}

double MD::initialise(vector_3d &r, vector_3d &v, double TEMPERATURE) {
  /* Initialise position matrix and velocity matrix from Cubic Centred Lattice
   */
  if (!options.compression_options.compression ||
      (options.compression_options.compression &&
       options.compression_options.compress_count == 0)) {
    /* Select the lattice formation */
    choose_lattice_formation(options.lattice, r);

    /* Generates Maxwell-Boltzmann distribution */
    mb_distribution(v, TEMPERATURE);
  }

  /* Calculate the average velocities */
  double mean_vx = std::accumulate(v.x.begin(), v.x.end(), 0.0) / options.N;
  double mean_vy = std::accumulate(v.y.begin(), v.y.end(), 0.0) / options.N;
  double mean_vz = std::accumulate(v.z.begin(), v.z.end(), 0.0) / options.N;
  /* Conserve the momentum of the fluid by subsracting the average velocities
     using a lambda expression */
  std::for_each(v.x.begin(), v.x.end(), [mean_vx](double &d) { d -= mean_vx; });
  std::for_each(v.y.begin(), v.y.end(), [mean_vy](double &d) { d -= mean_vy; });
  std::for_each(v.z.begin(), v.z.end(), [mean_vz](double &d) { d -= mean_vz; });

  size_t i;
  /* Temperature calculation, statistically */
  double KE = 0;
  for (i = 0; i < options.N; ++i) {
    KE += 0.5 * (v.x[i] * v.x[i] + v.y[i] * v.y[i] + v.z[i] * v.z[i]);
  }
  options.temperature = KE / (1.5 * options.N);
  options.scale_v =
      sqrt(TEMPERATURE / options.temperature); /* scaling factor */

  /* Velocity scaling */
  for (i = 0; i < options.N; ++i) {
    v.x[i] *= options.scale_v;
    v.y[i] *= options.scale_v;
    v.z[i] *= options.scale_v;
  }

  if (options.io_options.msd) {
    /* A copy of the r vectors where the BC will not be applied */
    MSD_r.x = r.x;
    MSD_r.y = r.y;
    MSD_r.z = r.z;

    /* MSD initialisation, storing first positions of particles */
    MSD.x = r.x;
    MSD.y = r.y;
    MSD.z = r.z;
  }

  if (options.io_options.vaf) {
    /* VAF initialisation, storing first velocities of particles */
    Cv.x = v.x;
    Cv.y = v.y;
    Cv.z = v.z;
  }

  return KE;
}

void MD::mb_distribution(vector_3d &v, double TEMPERATURE) {
  double kb = 1.0;
  double m = 1.0;

  double var = sqrt(TEMPERATURE * kb / m);
  double mean = 0;

  /* Use current time as seed for random generator */
  std::srand(std::time(nullptr));
  int random_variable = std::rand();
  if (options.test_options.is_testing)
    random_variable = 666; /* Fixing it for testing */

  std::default_random_engine generator;
  generator.seed(random_variable);

  std::normal_distribution<double> g_x(mean, var);
  std::normal_distribution<double> g_y(mean, var);
  std::normal_distribution<double> g_z(mean, var);

  for (size_t i = 0; i < options.N; ++i) {
    v.x.push_back(g_x(generator));
    v.y.push_back(g_y(generator));
    v.z.push_back(g_z(generator));
  }
}

double MD::verlet_algorithm(vector_3d &r, vector_3d &v, vector_3d &f,
                            bool sample_msd = true) {
  size_t i;
  double KE = 0;

  for (i = 0; i < options.N; ++i) {
    /*************************************************************************/
    /* Step velocities forward in time */
    v.x[i] = v.x[i] * options.scale_v + f.x[i] * options.dt;
    v.y[i] = v.y[i] * options.scale_v + f.y[i] * options.dt;
    v.z[i] = v.z[i] * options.scale_v + f.z[i] * options.dt;

    /* Step positions forward in time */
    r.x[i] = r.x[i] + v.x[i] * options.dt;
    r.y[i] = r.y[i] + v.y[i] * options.dt;
    r.z[i] = r.z[i] + v.z[i] * options.dt;

    if (sample_msd) {
      /* MSD stepping */
      MSD_r.x[i] += v.x[i] * options.dt;
      MSD_r.y[i] += v.y[i] * options.dt;
      MSD_r.z[i] += v.z[i] * options.dt;
    }
    /**************************************************************************/

    /* Kinetic Energy Calculation */
    KE += 0.5 * (v.x[i] * v.x[i] + v.y[i] * v.y[i] + v.z[i] * v.z[i]);

    /* Apply periodic boundary conditions to ensure particles remain
       inside the box */
    // todo: make boundary conditions routines
    if (r.x[i] > options.Lx) r.x[i] = r.x[i] - options.Lx;
    if (r.x[i] < 0.0) r.x[i] = r.x[i] + options.Lx;
    if (r.y[i] > options.Ly) r.y[i] = r.y[i] - options.Ly;
    if (r.y[i] < 0.0) r.y[i] = r.y[i] + options.Ly;
    if (r.z[i] > options.Lz) r.z[i] = r.z[i] - options.Lz;
    if (r.z[i] < 0.0) r.z[i] = r.z[i] + options.Lz;
  }

  return KE;
}

void MD::velocity_autocorrelation_function(vector_3d &Cv, vector_3d &v) {
  double cr_temp = 0; /* resets the sum every time step */
  double m = 1.0;     /* particle mass */
  size_t i;
  /* The peak of the VAF is located at 3kb*T/m */
  double norm = 3 * options.target_temperature / m;
  for (i = 0; i < options.N; i++) {
    cr_temp += (Cv.x[i] * v.x[i] + Cv.y[i] * v.y[i] + Cv.z[i] * v.z[i]);
  }
  Cr.push_back((cr_temp / options.N) / norm);
}

void MD::radial_distribution_function(double &rho, double &cut_off,
                                      size_t &bins, size_t &particles,
                                      std::ofstream &fstream) {
  double R = 0;
  double norm = 1;
  /* Exclude the self particle interaction from the density */
  double cor_rho = rho * (particles - 1) / particles;
  double dr = cut_off / bins;

  fstream << "# particles (N): " << particles << " steps: " << options.steps
          << " rho: " << rho << " bins: " << bins
          << " cut_off (rg): " << cut_off << " dr: " << dr << std::endl;
  fstream << "# Radius (r)" << '\t' << "Normalised" << '\t' << "Unormalised"
          << std::endl;

  for (size_t i = 1; i < bins; ++i) {
    R = cut_off * i / bins;
    /* Volume between 2 spheres, accounting for double counting
     hence the 2/3*pi*((R+dr)**3 - R**3)
     Accounting for the rdf_wait time steps */
    norm = cor_rho * (2.0 / 3.0 * PI * particles *
                      (options.steps - options.rdf_options.rdf_wait) *
                      (pow((R + (dr / 2.0)), 3) - pow((R - (dr / 2.0)), 3)));

    fstream << R << '\t' << rdf[i] / norm << '\t' << rdf[i] << std::endl;
  }
}

void MD::mean_square_displacement(vector_3d &MSD, vector_3d &MSD_r) {
  double msd_temp = 0;

  for (size_t i = 0; i < options.N; ++i) {
    msd_temp +=
        (pow((MSD_r.x[i] - MSD.x[i]), 2) + pow((MSD_r.y[i] - MSD.y[i]), 2) +
         pow((MSD_r.z[i] - MSD.z[i]), 2));
  }
  msd.push_back(msd_temp / options.N);
}

void MD::structure_factor(vector_3d &r) {
  double s = pow((options.N / options.density), (1.0 / 3.0));
  double fkx1 = 2.0 * PI / (s / (2.0 * options.Nx));
  double fky1 = 2.0 * PI / (s / (2.0 * options.Ny));
  double fkz1 = 2.0 * PI / (s / (2.0 * options.Nz));
  double sfcosx = 0, sfcosy = 0, sfcosz = 0;
  double sfsinx = 0, sfsiny = 0, sfsinz = 0;

  /* Try to calculate the structure factor at once for all axis
     if the particles per axis are equal. Else simply do them individually */

  if (options.Nx == options.Ny && options.Nx == options.Nz) {
    for (size_t i = 0; i < r.x.size(); ++i) {
      sfcosx += cos(fkx1 * r.x[i]);
      sfsinx += sin(fkx1 * r.x[i]);
      sfcosy += cos(fky1 * r.y[i]);
      sfsiny += sin(fky1 * r.y[i]);
      sfcosz += cos(fkz1 * r.z[i]);
      sfsinz += sin(fkz1 * r.z[i]);
    }
  } else {
    for (const auto &i : r.x) {
      sfcosx += cos(fkx1 * i);
      sfsinx += sin(fkx1 * i);
    }
    for (const auto &i : r.y) {
      sfcosy += cos(fky1 * i);
      sfsiny += sin(fky1 * i);
    }
    for (const auto &i : r.z) {
      sfcosz += cos(fkz1 * i);
      sfsinz += sin(fkz1 * i);
    }
  }

  double kx = sqrt(pow(sfcosx / options.N, 2) + pow(sfsinx / options.N, 2));
  double ky = sqrt(pow(sfcosy / options.N, 2) + pow(sfsiny / options.N, 2));
  double kz = sqrt(pow(sfcosz / options.N, 2) + pow(sfsinz / options.N, 2));

  sf.x.push_back(kx);
  sf.y.push_back(ky);
  sf.z.push_back(kz);
}

std::tuple<double, double> MD::calculate_forces(size_t &step_idx,
                                                pair_potential_type force) {
  /* Resetting forces */
  std::fill(f.x.begin(), f.x.end(), 0);
  std::fill(f.y.begin(), f.y.end(), 0);
  std::fill(f.z.begin(), f.z.end(), 0);

  /* Reseting <Potential> U to 0 */
  double U = 0;  /* Potential Energy */
  double PC = 0; /* Configurational Pressure */

  size_t i, j, igr;
  for (i = 0; i < options.N - 1; ++i) {
    for (j = i + 1; j < options.N; ++j) {
      /* distance between particle i and j */
      double x = r.x[i] - r.x[j]; /* Separation distance */
      double y = r.y[i] - r.y[j]; /* between particles i and j */
      double z = r.z[i] - r.z[j]; /* in Cartesian */

      /* Get the shortest image of the two particles
         if the particles are near the periodic boundary,
         this image is their reflection. */
      if (x > (0.5 * options.L)) x = x - options.L;
      if (x < (-0.5 * options.L)) x = x + options.L;
      if (y > (0.5 * options.L)) y = y - options.L;
      if (y < (-0.5 * options.L)) y = y + options.L;
      if (z > (0.5 * options.L)) z = z - options.L;
      if (z < (-0.5 * options.L)) z = z + options.L;

      /* Pair potential radius */
      double radius = sqrt((x * x) + (y * y) + (z * z));

      /* Force loop */
      if (radius < options.cut_off) {
        /* Allows the user to choose different pair potentials */
        auto [ff, temp_u] = force(radius, options.power, options.a_cst);

        /* Average potential energy */
        if (options.io_options.energies) U += temp_u;

        /* Configurational pressure */
        if (options.io_options.pressure) PC += radius * ff;

        /* Canceling the ij and ji pairs
           Taking the lower triangular matrix */
        f.x[i] += x * ff / radius;
        f.x[j] -= x * ff / radius;
        f.y[i] += y * ff / radius;
        f.y[j] -= y * ff / radius;
        f.z[i] += z * ff / radius;
        f.z[j] -= z * ff / radius;

        /* Radial Distribution
           measured with a delay, since the system requires a few thousand
           time-steps to reach equilibrium */
        if (options.io_options.rdf) {
          if (step_idx > options.rdf_options.rdf_wait) {
            igr =
                round(options.rdf_options.rdf_bins * radius / options.cut_off);
            rdf[igr] += 1;
          }
        }
      }
    }
  }

  return std::make_tuple(U, PC);
}

void MD::simulation(std::string simulation_name, double DENSITY,
                    double TEMPERATURE, double POWER, double A_CST,
                    std::string pp_type) {
  /* NOTE: this is a legacy routine and it will be removed in the future
     Initialise the variables with the input parameters
     Name the simulation. This will be used as a prefix in the files */
  options.io_options.simulation_name = simulation_name;
  options.density = DENSITY;
  options.target_temperature = TEMPERATURE;
  options.power = POWER;
  options.a_cst = A_CST;
  options.potential_type = pp_type;

  simulation();
}

void MD::simulation() {
  /* Preallocate storage space */
  set_vector_sizes();

  std::cout << "***************************\n"
               "** MD simulation started **\n"
               "***************************\n"
            << std::endl;
  /* Sets the unneeded variables (A and/or n) to NAN depending on the pp-type */
  std::cout << set_simulation_params(options.density,
                                     options.target_temperature, options.power,
                                     options.a_cst, options.potential_type)
            << std::endl;

  /* Gets the pair potential for the simulation based on a map of the
     initials of the pair potential and the pair potential itself. */
  pair_potential_type force = get_force_func(options.potential_type);

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  /* Initialise the simulation, lattice params and much more */
  options.kinetic_energy = initialise(r, v, options.target_temperature);
  size_t step_idx;
  for (step_idx = 0; step_idx < options.steps; ++step_idx) {
    /* Forces loop */
    /* Resetting forces to zero */
    std::fill(f.x.begin(), f.x.end(), 0);
    std::fill(f.y.begin(), f.y.end(), 0);
    std::fill(f.z.begin(), f.z.end(), 0);

    /* Reseting <Potential> U to 0 */
    double U = 0;  /* Potential Energy */
    double PC = 0; /* Configurational Pressure */

    size_t i, j, igr;
    for (i = 0; i < options.N - 1; ++i) {
      for (j = i + 1; j < options.N; ++j) {
        /* distance between particle i and j */
        double x = r.x[i] - r.x[j]; /* Separation distance */
        double y = r.y[i] - r.y[j]; /* between particles i and j */
        double z = r.z[i] - r.z[j]; /* in Cartesian */

        /* Get the shortest image of the two particles
           if the particles are near the periodic boundary,
           this image is their reflection. */
        if (x > (0.5 * options.Lx)) x -= options.Lx;
        if (x < (-0.5 * options.Lx)) x += options.Lx;
        if (y > (0.5 * options.Ly)) y -= options.Ly;
        if (y < (-0.5 * options.Ly)) y += options.Ly;
        if (z > (0.5 * options.Lz)) z -= options.Lz;
        if (z < (-0.5 * options.Lz)) z += options.Lz;

        /* Pair potential radius */
        double radius = sqrt((x * x) + (y * y) + (z * z));

        /* Force loop */
        if (radius < options.cut_off) {
          /* Allows the user to choose different pair potentials */
          auto [ff, temp_u] = force(radius, options.power, options.a_cst);

          /* Average potential energy */
          U += temp_u;

          /* Configurational pressure */
          PC += radius * ff;

          /* Canceling the ij and ji pairs
             Taking the lower triangular matrix */
          f.x[i] += x * ff / radius;
          f.x[j] -= x * ff / radius;
          f.y[i] += y * ff / radius;
          f.y[j] -= y * ff / radius;
          f.z[i] += z * ff / radius;
          f.z[j] -= z * ff / radius;

          /* Radial Distribution
             measured with a delay, since the system requires a few thousand
             time-steps to reach equilibrium */
          if (step_idx > options.rdf_options.rdf_wait) {
            igr =
                round(options.rdf_options.rdf_bins * radius / options.cut_off);
            rdf[igr] += 1;
          }
        }
      }
    }

    /* Isothermal Calibration */
    /* using T & KE from prev timestep */
    options.scale_v = sqrt(options.target_temperature / options.temperature);

    options.kinetic_energy = verlet_algorithm(r, v, f, true);

    if (options.io_options.msd) mean_square_displacement(MSD, MSD_r);

    if (options.io_options.vaf) velocity_autocorrelation_function(Cv, v);

    /* Calculate the structure factor k-vectors */
    if (options.io_options.sf) structure_factor(r);

    /* Average Temperature */
    options.temperature = options.kinetic_energy / (1.5 * options.N);
    temperature.push_back(options.temperature);

    /* Density */
    density.push_back(options.density);

    if (options.io_options.pressure) {
      /* Average Configurational Pressure Pc */
      pc.push_back(PC / (3 * options.volume));

      /* Kinetic Pressure */
      pk.push_back(options.density * options.temperature);
    }

    if (options.io_options.energies) {
      /* Average Potential Energy per particle */
      u_en.push_back(U / options.N);

      /* Average Kinetic Energy */
      k_en.push_back(options.kinetic_energy / options.N);
    }

    /* Save positions for visualisation with Python */
    if (options.io_options.visualise) {
      /* Reserve memory for the position vectors */
      (*pos_x)[step_idx].reserve(options.N);
      (*pos_y)[step_idx].reserve(options.N);
      (*pos_z)[step_idx].reserve(options.N);

      /* Populate the vectors with the current positions */
      (*pos_x)[step_idx] = r.x;
      (*pos_y)[step_idx] = r.y;
      (*pos_z)[step_idx] = r.z;
    }
  }
  /****************************************************************************/
  /* simulation Ends HERE */

  /* Write to files */
  file_output();

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

  /* Close file streams, makes simulation reusable in loops */
  reset_values();
}

void MD::reset_values(bool force_reset) {
  /*
    bc the stream might be closed and the user might then call reset_values
    which will throw an exception
  */

  /* Do not close streams and do not clear position and velocity vectors
     in the case where the fluid is being compressed */
  if (!options.compression_options.compression || force_reset) {
    /* Clear values, size, but reserve capacity */
    r.x.clear();
    r.y.clear();
    r.z.clear();
    v.x.clear();
    v.y.clear();
    v.z.clear();
    options.compression_options.compress_count = 0;
  }
  /* Reset the MSD initial vectors */
  MSD_r.x.clear();
  MSD_r.y.clear();
  MSD_r.z.clear();
  /* Clear monitored quantities */
  density.clear();
  temperature.clear();
  u_en.clear();
  k_en.clear();
  pc.clear();
  pk.clear();
  msd.clear();
  Cr.clear();
  rdf.resize(options.rdf_options.rdf_bins + 1, 0); /* gr with Index igr */
}

std::string MD::set_simulation_params(double &rho, double &T, double &power,
                                      double &a, std::string &pp_type) {
  std::string params =
      "Fluid parameters: rho: " + stat_file::convert_to_string(rho, 4) +
      " T: " + stat_file::convert_to_string(T, 4);

  params = "Lattice: " + options.lattice + "\n" + params;

  if (pp_type == "GaussianCoreModel") {
    params = "Potential: GaussianCoreModel, " + params;
    options.power =
        NAN; /* Set the variable to NAN to be ignore by the logger */
    options.a_cst =
        NAN; /* Set the variable to NAN to be ignore by the logger */
  }

  else if (pp_type == "LennardJones") {
    params = "Potential: LennardJones, " + params;
    options.power =
        NAN; /* Set the variable to NAN to be ignore by the logger */
    options.a_cst =
        NAN; /* Set the variable to NAN to be ignore by the logger */
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

void MD::enable_testing(bool is_testing) {
  options.test_options.is_testing = is_testing;
}

void MD::set_vector_sizes() {
  /* For efficiency, memory in the containers is reserved before use */
  /* Positions */
  r.x.reserve(options.N);
  r.y.reserve(options.N);
  r.z.reserve(options.N);
  /* Velocities */
  v.x.reserve(options.N);
  v.y.reserve(options.N);
  v.z.reserve(options.N);
  /* RDF */
  if (options.io_options.rdf)
    rdf.resize(options.rdf_options.rdf_bins + 1, 0); /* gr with Index igr */
  /* Forces/Acceleration */
  f.x.resize(options.N, 0);
  f.y.resize(options.N, 0);
  f.z.resize(options.N, 0);
  /* Structure factor k-arrays */
  if (options.io_options.sf) {
    sf.x.reserve(options.N);
    sf.y.reserve(options.N);
    sf.z.reserve(options.N);
  }
  /* Observed Quantities */
  if (options.io_options.vaf)
    Cr.reserve(options.steps); /* Velocity Autocorrelation Function */
  if (options.io_options.msd)
    msd.reserve(options.steps); /* Mean Square Displacement */
  if (options.io_options.energies) {
    u_en.reserve(options.steps); /* Average Potential Energy */
    k_en.reserve(options.steps); /* Average Kinetic Energy */
  }
  if (options.io_options.pressure) {
    pc.reserve(options.steps); /* Configurational Pressure */
    pk.reserve(options.steps); /* Kinetic Pressure */
  }
  temperature.reserve(options.steps);
}

void MD::save_visualisation_arrays() {
  /* Write the arrays as jagged,(hence transposed), this creates rows=STEPS */
  /* and columns=PARTICLES */
  std::ofstream out_x(
      logger.file_naming(options.io_options.dir + "/" +
                             options.io_options.simulation_name + "x_data",
                         options.steps, options.N, options.density,
                         options.target_temperature, options.power,
                         options.a_cst),
      std::ofstream::trunc | std::ofstream::out);
  logger.write_file(*pos_x, out_x, "");
  out_x.close();

  std::ofstream out_y(
      logger.file_naming(options.io_options.dir + "/" +
                             options.io_options.simulation_name + "y_data",
                         options.steps, options.N, options.density,
                         options.target_temperature, options.power,
                         options.a_cst),
      std::ofstream::trunc | std::ofstream::out);
  logger.write_file(*pos_y, out_y, "");
  out_y.close();

  std::ofstream out_z(
      logger.file_naming(options.io_options.dir + "/" +
                             options.io_options.simulation_name + "z_data",
                         options.steps, options.N, options.density,
                         options.target_temperature, options.power,
                         options.a_cst),
      std::ofstream::trunc | std::ofstream::out);
  logger.write_file(*pos_z, out_z, "");
  out_z.close();
}

void MD::file_output() {
  /* Generating the filenames for the output */

  /* Always open the data file stream since the temperature and density */
  /* are always written */
  logger.file_names.push_back(
      options.io_options.dir +
      logger.file_naming("/" + options.io_options.simulation_name + "Data",
                         options.steps, options.N, options.density,
                         options.target_temperature, options.power,
                         options.a_cst));

  if (options.io_options.position) {
    logger.file_names.push_back(
        options.io_options.dir +
        logger.file_naming(
            "/" + options.io_options.simulation_name + "Positions_Velocities",
            options.steps, options.N, options.density,
            options.target_temperature, options.power, options.a_cst));
  }

  if (options.io_options.rdf) {
    logger.file_names.push_back(
        options.io_options.dir +
        logger.file_naming("/" + options.io_options.simulation_name + "RDF",
                           options.steps, options.N, options.density,
                           options.target_temperature, options.power,
                           options.a_cst));
  }

  /* Save particle positions to files */
  if (options.io_options.visualise) save_visualisation_arrays();

  /* Create a map with all the streams */
  std::vector<std::ofstream> all_streams = logger.open_files(logger.file_names);
  std::map<std::string, std::ofstream> streams;
  streams["data"] = std::move(all_streams[0]);

  for (size_t stream = 1; stream < all_streams.size(); ++stream) {
    if (options.io_options.position && !streams.count("position"))
      streams["position"] = std::move(all_streams[stream]);
    else if (options.io_options.rdf && !streams.count("rdf"))
      streams["rdf"] = std::move(all_streams[stream]);
    else
      std::cerr << "Unrecognised stream in all_streams" << std::endl;
  }

  /* generate the correct header depending on io_options */
  std::vector<std::vector<double>> output_quantities = {density, temperature};
  std::string header = "# step\trho\tT";

  if (options.io_options.energies) {
    output_quantities.push_back(u_en);
    output_quantities.push_back(k_en);
    header += "\tU\tK";
  }
  if (options.io_options.pressure) {
    output_quantities.push_back(pc);
    output_quantities.push_back(pk);
    header += "\tPc\tPk";
  }
  if (options.io_options.msd) {
    output_quantities.push_back(msd);
    header += "\tMSD";
  }
  if (options.io_options.vaf) {
    output_quantities.push_back(Cr);
    header += "\tVAF";
  }
  if (options.io_options.sf) {
    output_quantities.push_back(sf.x);
    output_quantities.push_back(sf.y);
    output_quantities.push_back(sf.z);
    header += "\tSFx\tSFy\tSFz";
  }
  logger.write_data_file(streams["data"], header, output_quantities);

  /* Saving Last Position */
  if (options.io_options.position) {
    logger.write_data_file(streams["position"],
                           "# particle\tx\ty\tz\tvx\tvy\tvz\tax\tay\taz",
                           {r.x, r.y, r.z, v.x, v.y, v.z, f.x, f.y, f.z});
  }

  if (options.io_options.rdf) {
    radial_distribution_function(options.density, options.cut_off,
                                 options.rdf_options.rdf_bins, options.N,
                                 streams["rdf"]);
  }

  for (auto &[key, val] : streams) val.close();
}