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
#include <fstream>
#include <vector>

class MD {
 protected:
  std::vector<double> rx, ry, rz;        // Position Arrays
  std::vector<double> vx, vy, vz;        // Velocity Arrays
  std::vector<double> fx, fy, fz;        // Force arrays
  std::vector<double> Cvx, Cvy, Cvz;     // VAF arrays
  std::vector<double> rrx, rry, rrz;     // used in MSD calculation
  std::vector<double> MSDx, MSDy, MSDz;  // MSD arrays
  // Statistical quantity vectors, VAF, MSD, Energies and pressures
  std::vector<double> Cr, msd, u_en, k_en, pc, pk, temperature, density;

  size_t Nx, Ny, Nz, nhist, rdf_wait;  // Particles in x, y, z
  size_t N, _STEP_INDEX, STEPS;  // Total particles, step index, maximum steps
  double _T0;                    // Target/ Thermostat temperature
  double dt;                     // time step
  double x, y, z;                // distance between particle i and j
  double r;                      // distance in polar
  double _rho;                   // density
  double scale;                  // box scaling parameter
  double KE = 0.0;               // Kinetic Energy
  double T;                      // Temperature
  double L;                      // Length of the box after scaling
  double Vol;                    // Volume
  double cut_off = 3.0;          // simulation runs only within cutoff
  double U = 0;                  // Potential Energy
  double PC = 0;                 // Configurational Pressure
  double PK;                     // Kinetic Pressure
  double scale_v;                // velocity scaling

  // Visualisation vectors initialised in constructor
  std::vector<std::vector<double>> *pos_x;
  std::vector<std::vector<double>> *pos_y;
  std::vector<std::vector<double>> *pos_z;

  // Compression variables
  bool compress = false;
  size_t c_counter = 0;  // counts the number compression that have occurred

  // HISTOGRAM VARIABLES
  int igr;                 // Index of Hist
  double rg;               // cut off radius
  double dr;               // bin increment in terms of radius units
  std::vector<double> gr;  // RDF vector container

 private:
  long double PI;
  std::string _FILE_EXT;  // output file extension
  /* Variables for storing inside the object the file ID */
  std::string full_exe_dir, top_exe_dir;
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str,
      _n_to_str, _A_to_str;
  std::string rdf, data, pos;
  std::string _dir, _FILE_ID;
  std::ofstream RDF, DATA, POS;

 public:
  bool VISUALISE;
  MD(std::string DIRECTORY, size_t run_number);
  MD(std::string DIRECTORY, size_t run_number, bool COMPRESS_FLAG);
  MD(std::string DIRECTORY, size_t run_number, bool COMPRESS_FLAG,
     size_t rdf_bins, size_t particles_per_axis, bool track_particles,
     size_t collect_rdf_after);

  ~MD();

  void Simulation(double DENSITY, double TEMPERATURE, int POWER, double A_CST);
  void get_phases(double DENSITY, double FINAL_DENSITY, double DENSITY_INC,
                  double TEMPERATURE, int POWER, double A_CST);
  void reset_values();

 protected:
  void initialise(std::vector<double> &x, std::vector<double> &y,
                  std::vector<double> &z, std::vector<double> &vx,
                  std::vector<double> &vy, std::vector<double> &vz,
                  double TEMPERATURE);
  void mb_distribution(double TEMPERATURE);
  void verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry,
                        std::vector<double> &rz, std::vector<double> &vx,
                        std::vector<double> &vy, std::vector<double> &vz,
                        std::vector<double> &rrx, std::vector<double> &rry,
                        std::vector<double> &rrz);
  void velocity_autocorrelation_function(std::vector<double> &Cvx,
                                         std::vector<double> &Cvy,
                                         std::vector<double> &Cvz);
  void radial_distribution_function();
  void mean_square_displacement(std::vector<double> &MSDx,
                                std::vector<double> &MSDy,
                                std::vector<double> &MSDz);

  void open_files();
  std::string file_naming(std::string prefix, double DENSITY,
                          double TEMPERATURE, int POWER, double A_cst);
  void write_to_files();
  void show_run(size_t step_size_show);
  void time_stamp(std::ofstream &, std::string variables);
  std::string convert_to_string(const double &x, const int &precision);
};