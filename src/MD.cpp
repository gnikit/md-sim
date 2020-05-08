#include "MD.h"
// todo: add logger https://github.com/gabime/spdlog
// TODO: scale the box by Lx, Ly, Lz in a tensor form

MD::MD() {}

MD::MD(options_type &input_options) {
  /* Perform a shallow copy of the input_options to options */
  options = input_options;

  /* Test whether the input directory exists */
  if (!input_options.io_options.dir.empty()) {
    try {
      if (!fs::exists(input_options.io_options.dir)) {
        throw
          "input out_directory in MD constructor does not exist.\n"
          "Use a valid directory for output files to be saved";
      }

    } catch (const char *msg) {
      std::cerr << "Error: " << msg << std::endl;
      exit(-1);
    }
  }

  std::cout << "Output directory set to: " << options.io_options.dir
            << std::endl;

  /* Print type of simulation */
  std::cout << "Simulation type: " << options.simulation_type << std::endl;

  /* Print simulation name if any */
  std::cout << "Simulation name: " << options.io_options.simulation_name
            << std::endl;

  /* Are we saving all the positions for the fluid */
  std::cout << "Particle visualisation: " << options.io_options.visualise
            << std::endl;

  /* Print stepping algorithm */
  std::cout << "Iterative algorithm: " << options.iterative_method << std::endl;

  /* Print number of iterations */
  std::cout << "Number of steps: " << options.steps << std::endl;

  /* Print particles and lattice */
  std::cout << "Initial lattice: " << options.lattice << std::endl;

  /* Check particles have been supplied in the correct form */
  try {
    if (input_options.particles.empty()) {
      throw "The supplied particles vector is empty";
    } else if (input_options.particles.size() < 3) {
      throw "The supplied particles vector must be of size 3";
    } else if (std::find(input_options.particles.begin(),
                         input_options.particles.end(),
                         0) != input_options.particles.end()) {
      throw
            "The supplied particles vector contains a 0\n"
            "particles cannot be 0 in x, y or z";
    }

  } catch (const char *msg) {
    std::cerr << "Error: " << msg << std::endl;
    exit(-1);
  }

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
  /* Print the pair potential */
  std::cout << "Pair potential: " << options.potential_type << std::endl;

  if (input_options.density < 0) {
    std::cerr << "Error: Negative density supplied" << std::endl;
    exit(-1);
  }

  if (input_options.target_temperature < 0) {
    std::cerr << "Error: Negative temperature supplied" << std::endl;
    exit(-1);
  }

  /* Initialise scaling variables */
  options.dt = 0.005 / sqrt(options.target_temperature);  // todo: add to schema
  /* Box length scaling */
  options.L = pow((options.N / options.density), 1.0 / 3.0);
  options.Lx = options.Ly = options.Lz = options.L;  // todo: questionable!
  options.volume = options.N / options.density;

  /* cut_off definition */
  if (input_options.cut_off > 0) {
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

  /* Set boundary conditions */
  std::cout << "Boundary conditions: " << options.bcs << std::endl;

  /* Accuracy of RDF */
  std::cout << "RDF bins: " << options.rdf_options.rdf_bins << std::endl;

  /* Ensuring the number of steps is greater than the rdf equilibration period
   */
  try {
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
  /* Print testing options */
  std::cout << "Testing: " << options.test_options.is_testing << std::endl;

  /* Pass the options reference back to input_options to update the variable.
   * Routines in the phase_transition class depend on this line */
  input_options = options;

  /* Visualisation vectors on the heap*/
  pos = new std::vector<std::vector<double>>(4);
  for (size_t i = 0; i < (*pos).size(); ++i) (*pos)[i].reserve(options.N);
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
  r.reserve(options.N);
  /* Velocities */
  v.reserve(options.N);
  /* RDF */
  rdf.resize(options.rdf_options.rdf_bins + 1, 0); /* gr with Index igr */
  /* Forces/Acceleration */
  f.resize(options.N, 0);
  /* Structure factor k-arrays */
  sf.reserve(options.N);
  /* Observed Quantities */
  Cr.reserve(options.steps);   /* Velocity Autocorrelation Function */
  msd.reserve(options.steps);  /* Mean Square Displacement */
  u_en.reserve(options.steps); /* Average Potential Energy */
  k_en.reserve(options.steps); /* Average Kinetic Energy */
  pc.reserve(options.steps);   /* Configurational Pressure */
  pk.reserve(options.steps);   /* Kinetic Pressure */
  temperature.reserve(options.steps);

  /* Visualisation vectors on the heap*/
  pos = new std::vector<std::vector<double>>(4);
  for (size_t i = 0; i < (*pos).size(); ++i) (*pos)[i].reserve(options.N);
  options.test_options.is_testing = false;
}

/* Delegating constructors with reduced number of arguments
 https://en.wikipedia.org/wiki/C++11#Object_construction_improvement
 Convinient constructor to use for simple cases
 */

MD::~MD() {
  /* Destroy the vectors allocated on the heap */
  /* There is no pointer to delete if there are no particles*/
  // todo: use vector of smart pointers
  if (options.N != 0) delete pos;
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

  auto begin = std::chrono::high_resolution_clock::now();

  /* Initialise the simulation, lattice params and much more */
  options.kinetic_energy = initialise(r, v, options.target_temperature);
  size_t step_idx;
  for (step_idx = 0; step_idx < options.steps; ++step_idx) {
    /* Isothermal Calibration */
    /* using T & KE from prev timestep */
    options.scale_v = sqrt(options.target_temperature / options.temperature);

    double U, PC;
    std::tie(options.kinetic_energy, U, PC) =
        stepping_algorithm(r, v, f, step_idx, options.io_options.msd);

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
      /* Pass the pointer to the positions 2D vector */
      (*pos)[0] = r.x;
      (*pos)[1] = r.y;
      (*pos)[2] = r.z;
      (*pos)[3] = v.magnitude(); /* Speed for scaling */
      save_visualisation_arrays(step_idx);
    }
  }
  /****************************************************************************/
  /* simulation Ends HERE */

  /* Calculate the statistics and store them in a vector */
  run_stats = calculate_run_stats();

  /* Write to files */
  file_output(logger);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> sim_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);

  std::cout << "CPU run time = " << sim_time.count() << "s" << std::endl;
  std::cout << "******************************\n"
               "** MD simulation terminated **\n"
               "******************************\n"
            << std::endl;

  /* Close file streams, makes simulation reusable in loops */
  reset_values();
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

void MD::reset_values(bool force_reset) {
  /*
    bc the stream might be closed and the user might then call reset_values
    which will throw an exception
  */

  /* Do not close streams and do not clear position and velocity vectors
     in the case where the fluid is being compressed */
  if (!options.compression_options.compression || force_reset) {
    /* Clear values, size, but reserve capacity */
    r.clear();
    v.clear();
    options.compression_options.compress_count = 0;
  }
  /* Reset the MSD initial vectors */
  MSD_r.clear();
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

void MD::enable_testing(bool is_testing) {
  options.test_options.is_testing = is_testing;
}

/********************************** PROTECTED *********************************/

/*************************** INITIALISATION METHODS ***************************/

double MD::initialise(vector_3d<double> &r, vector_3d<double> &v,
                      double TEMPERATURE) {
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

  /* A copy of the r vectors where the BC will not be applied */
  /* MSD initialisation, storing first positions of particles */
  if (options.io_options.msd) MSD_r = MSD = r;

  /* VAF initialisation, storing first velocities of particles */
  if (options.io_options.vaf) Cv = v;

  return KE;
}

void MD::choose_lattice_formation(std::string &lattice, vector_3d<double> &r) {
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
    r.resize(options.N);
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

void MD::mb_distribution(vector_3d<double> &v, double TEMPERATURE) {
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

/***************************** ITERATIVE METHODS ******************************/

std::tuple<double, double, double> MD::stepping_algorithm(vector_3d<double> &r,
                                                          vector_3d<double> &v,
                                                          vector_3d<double> &f,
                                                          size_t &step,
                                                          bool msd) {
  double U = 0;
  double KE = 0;
  double PC = 0;

  /*************************************************************************/
  if (options.iterative_method == "Verlet") {
    std::tie(KE, U, PC) = verlet(r, v, f, step);

  } else if (options.iterative_method == "VelocityVerlet") {
    std::tie(KE, U, PC) = velocity_verlet(r, v, f, step);

  } else {
    std::cerr << "Iterative method: " << options.iterative_method
              << " has not been implemented. Aborting!" << std::endl;
    abort();
  }

  /* MSD stepping uses operator overlads */
  if (msd) MSD_r += v * options.dt;
  /**************************************************************************/

  apply_boundary_conditions(r, v, f);

  return std::make_tuple(KE, U, PC);
}

std::tuple<double, double, double> MD::verlet(vector_3d<double> &r,
                                              vector_3d<double> &v,
                                              vector_3d<double> &f,
                                              size_t &step) {
  size_t i;
  double KE = 0;
  double dt = options.dt;
  pair_potential_type force = get_force_func(options.potential_type);

  /* Calculate particle interactions based on input pair potential */
  auto [U, PC] = calculate_forces(r, f, step, force);

  for (i = 0; i < options.N; ++i) {
    /* v^(n+1) = v^(n)*scaled_v + f^(n)*dt */
    v.x[i] = v.x[i] * options.scale_v + f.x[i] * dt;
    v.y[i] = v.y[i] * options.scale_v + f.y[i] * dt;
    v.z[i] = v.z[i] * options.scale_v + f.z[i] * dt;

    /* r^(n+1) = r^(n) + v^(n+1)*dt */
    r.x[i] += v.x[i] * dt;
    r.y[i] += v.y[i] * dt;
    r.z[i] += v.z[i] * dt;

    /* Kinetic Energy Calculation */
    KE += 0.5 * (v.x[i] * v.x[i] + v.y[i] * v.y[i] + v.z[i] * v.z[i]);
  }

  return std::make_tuple(KE, U, PC);
}

std::tuple<double, double, double> MD::velocity_verlet(vector_3d<double> &r,
                                                       vector_3d<double> &v,
                                                       vector_3d<double> &f,
                                                       size_t &step) {
  double KE = 0;
  double dt = options.dt;
  vector_3d<double> f_n = f;
  pair_potential_type force = get_force_func(options.potential_type);

  /* r^(n+1) = r^(n) + v^(n)*dt + f^(n)*(dt^2 * 0.5) */
  for (size_t i = 0; i < options.N; ++i) {
    r.x[i] += v.x[i] * dt + 0.5 * f_n.x[i] * dt * dt;
    r.y[i] += v.y[i] * dt + 0.5 * f_n.y[i] * dt * dt;
    r.z[i] += v.z[i] * dt + 0.5 * f_n.z[i] * dt * dt;
  }

  /* f^(n+1) = calculate forces */
  auto [U, PC] = calculate_forces(r, f, step, force);

  /* v^(n+1) = v^(n) + (f^(n)+f^(n+1))*(dt*0.5) */
  for (size_t i = 0; i < options.N; ++i) {
    v.x[i] = v.x[i] * options.scale_v + (f.x[i] + f.x[i]) * (dt * 0.5);
    v.y[i] = v.y[i] * options.scale_v + (f.y[i] + f.y[i]) * (dt * 0.5);
    v.z[i] = v.z[i] * options.scale_v + (f.z[i] + f.z[i]) * (dt * 0.5);

    KE += 0.5 * (v.x[i] * v.x[i] + v.y[i] * v.y[i] + v.z[i] * v.z[i]);
  }

  return std::make_tuple(KE, U, PC);
}

std::tuple<double, double> MD::calculate_forces(vector_3d<double> &x,
                                                vector_3d<double> &f,
                                                size_t &step_idx,
                                                pair_potential_type &p) {
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
        auto [ff, temp_u] = p(radius, options.power, options.a_cst);

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
        if (options.io_options.rdf && step_idx > options.rdf_options.rdf_wait) {
          igr = round(options.rdf_options.rdf_bins * radius / options.cut_off);
          rdf[igr] += 1;
        }
      }
    }
  }

  return std::make_tuple(U, PC);
}

/***************************** BOUNDARY CONDITIONS ****************************/

void MD::apply_boundary_conditions(vector_3d<double> &r, vector_3d<double> &v,
                                   vector_3d<double> &f) {
  if (options.bcs == "Periodic") {
    periodic_boundary_conditions(r);

  } else if (options.bcs == "Reflective") {
    reflective_boundary_conditions(r, v);

  } else {
    std::cerr << "Supplied boundary condition: " << options.bcs
              << " not implemented. Aboring!" << std::endl;
    exit(-1);
  }
}

void MD::periodic_boundary_conditions(vector_3d<double> &r) {
  for (size_t i = 0; i < options.N; ++i) {
    if (r.x[i] > options.Lx) r.x[i] -= options.Lx;
    if (r.x[i] < 0.0) r.x[i] += options.Lx;
    if (r.y[i] > options.Ly) r.y[i] -= options.Ly;
    if (r.y[i] < 0.0) r.y[i] += options.Ly;
    if (r.z[i] > options.Lz) r.z[i] -= options.Lz;
    if (r.z[i] < 0.0) r.z[i] += options.Lz;
  }
}

void MD::reflective_boundary_conditions(vector_3d<double> const &r,
                                        vector_3d<double> &v) {
  for (size_t i = 0; i < options.N; ++i) {
    /* Flip the velocity sign */
    if (r.x[i] > options.Lx || r.x[i] < 0.0) v.x[i] = -v.x[i];
    if (r.y[i] > options.Ly || r.y[i] < 0.0) v.y[i] = -v.y[i];
    if (r.z[i] > options.Lz || r.z[i] < 0.0) v.z[i] = -v.z[i];
  }
}

/*********************** STATISTICAL QUANTITIES METHODS ***********************/

void MD::velocity_autocorrelation_function(vector_3d<double> &Cv,
                                           vector_3d<double> &v) {
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

  fstream << "#bins:" << bins << ",cut_off (rg):" << cut_off << ",dr:" << dr
          << std::endl;
  fstream << "#Radius (r),Normalised,Unormalised" << std::endl;

  for (size_t i = 1; i < bins; ++i) {
    R = cut_off * i / bins;
    /* Volume between 2 spheres, accounting for double counting
     hence the 2/3*pi*((R+dr)**3 - R**3)
     Accounting for the rdf_wait time steps */
    norm = cor_rho * (2.0 / 3.0 * PI * particles *
                      (options.steps - options.rdf_options.rdf_wait) *
                      (pow((R + (dr / 2.0)), 3) - pow((R - (dr / 2.0)), 3)));

    fstream << R << ',' << rdf[i] / norm << ',' << rdf[i] << std::endl;
  }
}

void MD::mean_square_displacement(vector_3d<double> &MSD,
                                  vector_3d<double> &MSD_r) {
  double msd_temp = 0;

  for (size_t i = 0; i < options.N; ++i) {
    msd_temp +=
        (pow((MSD_r.x[i] - MSD.x[i]), 2) + pow((MSD_r.y[i] - MSD.y[i]), 2) +
         pow((MSD_r.z[i] - MSD.z[i]), 2));
  }
  msd.push_back(msd_temp / options.N);
}

void MD::structure_factor(vector_3d<double> &r) {
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

/****************************** LOGGING METHODS *******************************/

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

void MD::save_visualisation_arrays(size_t dump_no) {
  std::string fname = options.io_options.dir + "/" +
                      options.io_options.simulation_name + "xyz_data_" +
                      std::to_string(dump_no) + ".csv";
  std::ofstream out_xyz(fname, std::ofstream::trunc | std::ofstream::out);
  logger.write_file(*pos, out_xyz, "x-pos, y-pos, z-pos, speed");
  out_xyz.close();
}

void MD::file_output(stat_file &logger) {
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
  std::string header = "#step,rho,T";

  if (options.io_options.energies) {
    // todo: I should be passing the pointer to u_en rather than the values
    output_quantities.push_back(u_en);
    output_quantities.push_back(k_en);
    header += ",U,K";
  }
  if (options.io_options.pressure) {
    output_quantities.push_back(pc);
    output_quantities.push_back(pk);
    header += ",Pc,Pk";
  }
  if (options.io_options.msd) {
    output_quantities.push_back(msd);
    header += ",MSD";
  }
  if (options.io_options.vaf) {
    output_quantities.push_back(Cr);
    header += ",VAF";
  }
  if (options.io_options.sf) {
    output_quantities.push_back(sf.x);
    output_quantities.push_back(sf.y);
    output_quantities.push_back(sf.z);
    header += ",SFx,SFy,SFz";
  }
  logger.write_data_file(streams["data"], header, output_quantities);

  /* Saving Last Position */
  if (options.io_options.position) {
    logger.write_data_file(streams["position"],
                           "#particle,x,y,z,vx,vy,vz,ax,ay,az",
                           {r.x, r.y, r.z, v.x, v.y, v.z, f.x, f.y, f.z});
    // todo: for performance I think this should be 2D pointer of pointers
  }

  if (options.io_options.rdf) {
    radial_distribution_function(options.density, options.cut_off,
                                 options.rdf_options.rdf_bins, options.N,
                                 streams["rdf"]);
  }

  /* Close all the streams */
  for (auto &[key, val] : streams) val.close();
  /* Empty the file_names generated since a compression calls the MD
     constructor once, hence `logger` is created only once */
  logger.file_names.clear();
}

std::vector<double> MD::calculate_run_stats() {
  /* Pass density as the first entry */
  std::vector<double> stat{options.density};

  /* Mean Square Displacement */
  if (options.io_options.msd) {
    std::vector<double> temp = vector_stats(msd);
    stat.insert(stat.end(), temp.begin(), temp.end());
  }

  if (options.io_options.vaf) {
    /* Velocity Autocorrelation Function */
    std::vector<double> temp = vector_stats(Cr);
    stat.insert(stat.end(), temp.begin(), temp.end());
  }

  /* Structure Factor */
  if (options.io_options.sf) {
    std::vector<double> temp = vector_stats(sf.x);
    stat.insert(stat.end(), temp.begin(), temp.end());
    temp = vector_stats(sf.y);
    stat.insert(stat.end(), temp.begin(), temp.end());
    temp = vector_stats(sf.z);
    stat.insert(stat.end(), temp.begin(), temp.end());
  }

  /* Potential Energy */
  if (options.io_options.energies) {
    std::vector<double> temp = vector_stats(u_en);
    stat.insert(stat.end(), temp.begin(), temp.end());
  }

  /* Configurational Pressure */
  if (options.io_options.pressure) {
    std::vector<double> temp = vector_stats(pc);
    stat.insert(stat.end(), temp.begin(), temp.end());
  }

  return stat;
}

std::vector<double> MD::get_run_stats() { return run_stats; }

/******************************** MISCELLANEOUS *******************************/

void MD::set_vector_sizes() {
  /* For efficiency, memory in the containers is reserved before use */
  /* Positions */
  r.reserve(options.N);
  /* Velocities */
  v.reserve(options.N);
  /* RDF */
  if (options.io_options.rdf)
    rdf.resize(options.rdf_options.rdf_bins + 1, 0); /* gr with Index igr */
  /* Forces/Acceleration */
  f.resize(options.N, 0);

  /* Observed Quantities */
  if (options.io_options.sf) sf.reserve(options.N); /* Structure Factor */
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