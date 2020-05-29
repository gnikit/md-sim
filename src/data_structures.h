#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>

/* Load Intel math lib if available */
#if defined(__INTEL_COMPILER)
#include <mathimf.h> /* Intel Math library */
#define COMPILER "INTEL"
#else
#include <math.h>
#endif

#include "vector_3d.h"
#include "vector_arithmetic_operators.h"

using namespace std;

struct rdf_options_type {
  size_t rdf_bins = 100;
  size_t rdf_wait = 0;
};
struct compression_options_type {
  bool compression = false;
  double density_final = 0;
  double density_inc = 0;
  bool reverse_comp = false;
  size_t compress_count = 0;
};

struct boundary_conditions_options {
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
  size_t verbosity = 0;          /* how verbose the standard output is */
  bool msd = true;               /* mean square displacement output */
  bool rdf = true;               /* radial distribution function output */
  bool vaf = true;               /* velocity autocorrelation output */
  bool energies = true;          /* potential energy output */
  bool pressure = true;          /* configurational pressure output */
  bool position = true;          /* particles' last positions output */
  bool sf = true;                /* structure factor output */
  bool visualise = false;        /* save all positions, of all particles */
  bool compression_stats = true; /* create a separate log file with stats */
  string dir = ".";              /* file output directory */
  string simulation_name = "";   /* simulation prefix/ name */
};

struct options_type {
  size_t dimension = 3;               /* geometrical dimension */
  string simulation_type = "";        /* type of simulation e.g. NormalRun */
  string potential_type = "";         /* pair potential type */
  string lattice = "SC";              /* lattice formation */
  string iterative_method = "Verlet"; /* iterative algorithm of particles */
  string bcs = "Periodic";            /* boundary conditions for whole box */
  vector<size_t> particles;           /* number of particles in each axis */
  size_t steps = 2000;                /* number of total iterations */
  size_t Nx, Ny, Nz, N = 0;           /* Particles in the x, y, z and total */
  double Lx, Ly, Lz, L = 0;           /* Individual box lengths */
  double volume = 0;                  /* volume of the box */
  double dt = 0.005;                  /* timestep */
  bool normalise_dt_w_temp = true;    /* normalise the timestep with T0 */
  double density = 0.5;               /* density */
  double target_temperature = 1.0;    /* target/ Thermostat temperature */
  double temperature = 1.0;           /* simulation temperature */
  double power = 0;                   /* pair potential intensity */
  double a_cst = 0;                   /* generic softening parameter */
  double q = 2.0;                     /* Parameter q for BIP */
  double kinetic_energy = 0;          /* kinetic energy */
  double cut_off = 0;                 /* cut off radius of simulation */
  double scale_v = 0;                 /* velocity scaling */
  bool fix_box_lengths = false;       /* Whether or not the MD cell has fix L */

  io_options_type io_options;
  rdf_options_type rdf_options;
  compression_options_type compression_options;
  test_options_type test_options;
};

class point_3d {
 public:
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
  vector_3d<double> r;     /* Position Arrays */
  vector_3d<double> v;     /* Velocity Arrays */
  vector_3d<double> f;     /* Force arrays */
  vector_3d<double> Cv;    /* VAF arrays */
  vector_3d<double> MSD_r; /* used in MSD calculation */
  vector_3d<double> MSD;   /* MSD arrays */
  vector_3d<double> sf;    /* Structure factor k-arrays */

  vector<double> Cr;          /* Velocity Autocorrelation Function*/
  vector<double> msd;         /* Mean Square Displacement */
  vector<double> u_en;        /*Potential Energy*/
  vector<double> k_en;        /*Kinetic Energy*/
  vector<double> pc;          /*Configurational Pressure*/
  vector<double> pk;          /* Kinetic Pressure*/
  vector<double> temperature; /*Temperature*/
  vector<double> density;     /*Density*/
  vector<double> rdf;         /*Radial distribution function*/

  options_type options; /*Options variables*/
};