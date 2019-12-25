
#include <map>
#include <string>
#include <vector>

using namespace std;

struct rdf_options_type {
  size_t rdf_bins = 200;
  size_t rdf_wait = 0;
};
struct compression_options_type {
  bool compression = false;
  double density_final = 0;
  double density_inc = 0;
  bool reverse_comp = false;
  size_t compress_count = 0;
};

struct boundary_contions_options_type {
  /**
   * @brief An array of boundaries where the map holds the type
   * of the boundary, e.g. Reflective, HardWall,
   * while the 2d array holds the points for the 2D surface
   */
  vector<map<string, vector<vector<double>>>> bc;
};

struct test_options_type {
  bool is_testing = false; /* fixes the seed for regression testing */
};

struct io_options_type {
  bool msd = true;             /* mean square displacement output */
  bool rdf = true;             /* radial distribution function output */
  bool vaf = true;             /* velocity autocorrelation output */
  bool energies = true;        /* potential energy output */
  bool pressure = true;        /* configurational pressure output */
  bool position = true;        /* particles' last positions output */
  bool sf = true;              /* structure factor output */
  bool visualise = false;      /* save all positions, of all particles */
  string dir = ".";            /* file output directory */
  string simulation_name = ""; /* simulation prefix/ name */

  io_options_type() {}
  io_options_type& operator=(io_options_type const& rhs) {
    msd = rhs.msd;
    rdf = rhs.rdf;
    vaf = rhs.vaf;
    energies = rhs.energies;
    pressure = rhs.pressure;
    position = rhs.position;
    sf = rhs.sf;
    visualise = rhs.visualise;
    dir = rhs.dir;
    simulation_name = rhs.simulation_name;

    return *this;
  }
};

struct options_type {
  string simulation_type = "";       /* type of simulation e.g. NormalRun */
  string potential_type = "";        /* pair potential type */
  string lattice = "";               /* lattice formation */
  string iterative_method = "";      /* iterative algorithm of particles */
  vector<size_t> particles{1, 1, 1}; /* number of particles in each axis */
  double random_lattice_var = 0;     /* variance of dist for random lattice */
  size_t steps = 2000;               /* number of total iterations */
  size_t Nx, Ny, Nz, N = 0;          /* Particles in the x, y, z and total */
  double Lx, Ly, Lz, L = 0;          /* Individual box lengths */
  double volume = 0;                 /* volume of the box */
  double dt = 0.005;                 /* timestep */
  bool normalise_dt_w_temp = true;   /* normalise the timestep with T0 */
  double density = 0.5;              /* density */
  double target_temperature = 0;     /* target/ Thermostat temperature */
  double temperature = 0;            /* simulation temperature */
  double power = 0;                  /* pair potential intensity */
  double a_cst = 0;                  /* generic softening parameter */
  double kinetic_energy = 0;         /* kinetic energy */
  double cut_off = 0;                /* cut off radius of simulation */
  double scale_v = 0;                /* velocity scaling */

  io_options_type io_options;
  rdf_options_type rdf_options;
  compression_options_type compression_options;
  test_options_type test_options;

  /* Default constructor */
  options_type() {}

  /* Overloaded assignment operator */
  options_type& operator=(options_type const& rhs) {
    simulation_type = rhs.simulation_type;
    potential_type = rhs.potential_type;
    lattice = rhs.lattice;
    iterative_method = rhs.iterative_method;

    random_lattice_var = rhs.random_lattice_var;
    steps = rhs.steps;

    density = rhs.density;
    target_temperature = rhs.target_temperature;
    power = rhs.power;
    a_cst = rhs.a_cst;

    cut_off = rhs.cut_off;

    /* Passing RDF options */
    rdf_options.rdf_bins = rhs.rdf_options.rdf_bins;
    rdf_options.rdf_wait = rhs.rdf_options.rdf_wait;

    /* testing options */
    test_options.is_testing = rhs.test_options.is_testing;

    /* Passing the io options */
    io_options = rhs.io_options;

    return *this;
  }
};

struct vector_3d {
  vector<double> x;
  vector<double> y;
  vector<double> z;
};

struct point_3d {
  double x;
  double y;
  double z;

  point_3d(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  inline point_3d operator=(const point_3d& a) {
    x = a.x;
    y = a.y;
    z = a.z;
    return *this;
  }

  inline point_3d operator+(const point_3d& a) const {
    return point_3d(x + a.x, y + a.y, z + a.z);
  }

  inline point_3d operator-(const point_3d& a) const {
    return point_3d(x - a.x, y - a.y, z - a.z);
  };

  inline point_3d operator*(const point_3d& a) const {
    return point_3d(x * a.x, y * a.y, z * a.z);
  }

  inline point_3d operator/(const point_3d& a) const {
    return point_3d(x / a.x, y / a.y, z / a.z);
  }
};

struct data_type {
  vector_3d r;     /* Position Arrays */
  vector_3d v;     /* Velocity Arrays */
  vector_3d f;     /* Force arrays */
  vector_3d Cv;    /* VAF arrays */
  vector_3d MSD_r; /* used in MSD calculation */
  vector_3d MSD;   /* MSD arrays */
  vector_3d sf;    /* Structure factor k-arrays */

  vector<double> Cr;
  vector<double> msd;
  vector<double> u_en;
  vector<double> k_en;
  vector<double> pc;
  vector<double> pk;
  vector<double> temperature;
  vector<double> density;
  vector<double> rdf;

  options_type options;
};