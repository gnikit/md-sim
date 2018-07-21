//////////////////////////////////////////////////////////////////////
// Ioannis Nikiteas 13/7/2017                                       //
//                                                                  //
// BSc Dissertation:                                                //
//  Investigating the transition from Molecular Dynamics to         //
//	Smoothed Particle Hydrodynamics                                 //
//                                                                  //
//	University: Royal Holloway University of London                 //
//                                                                  //
//	A program meant to simulate a MD fluid with an only             //
//	repulsive pair-potential. Increasing the parameter A            //
//	creates a coarse-graining effect for the system allowing it     //
//	to transition to SPH                                            //
//                                                                  //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#pragma once
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#elif defined(__GNUC__)
#include <math.h>
#define COMPILER "G++"
#else
#include <math.h>
#define COMPILER "WHO CARES"
#endif

#include <chrono>  // CPU run-time
#include <cstdint>
#include <ctime>    // std::chrono
#include <fstream>  // file writing
#include <iomanip>  // setprecision
#include <iostream>
#include <iterator>
#include <sstream>  // stringstream
#include <vector>

class MD {
 protected:
  //typedef std::vector<double> std::vector<double>;	// short-hand notation, use with care

  std::vector<double> rx, ry, rz;        // Position Arrays
  std::vector<double> vx, vy, vz;        // Velocity Arrays
  std::vector<double> fx, fy, fz;        // Force arrays
  std::vector<double> Cvx, Cvy, Cvz;     // VAF arrays
  std::vector<double> rrx, rry, rrz;     // used in MSD calculation
  std::vector<double> MSDx, MSDy, MSDz;  // MSD arrays
  std::vector<double> Cr, msd, u_en, k_en,
      pc, pk, temperature, density;  // vectors used as buffers

  size_t Nx, Ny, Nz;              // Particles in x, y, z
  size_t N, _STEP_INDEX, _STEPS;  // Total particles, step index, maximum steps
  double _T0;                     // Target Temperature. Desired T for the system to operate
  double dt = 0.005;              // time step This applies: dt = 0.005/sqrt(T0)
  double x, y, z;                 // distance between particle i and j
  double r;                       // distance in polar
  double _rho;                    // density
  double scale;                   // box scaling parameter
  double KE = 0.0;                // Kinetic Energy
  double T;                       // Temperature
  double L;                       // Length of the box after scaling
  double Vol;                     // Volume
  double cut_off = 3.0;           // simulation runs only within cutoff
  double U = 0;                   // Potential Energy
  double PC = 0;                  // Configurational Pressure
  double PK;                      // Kinetic Pressure
  double scale_v;                 // velocity scaling

  // Quenching varibles
  bool compression_flag = false;
  size_t Q_counter = 0;  // counts the number qunchings that have occured

  // HISTOGRAM VARIABLES
  int igr;  // Index of Hist
  double rg;
  double dr;
  std::vector<double> gr;  // RDF vector container

 private:
  const long double PI = acos(-1.0);
  std::string _FILE_EXT;  // output file extension
  /* Variables for storing inside the object the file ID */
  std::string full_exe_dir, top_exe_dir;
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str, _n_to_str, _A_to_str;
  std::string HIST, data, pos;
  std::string _dir, _FILE_ID;
  std::ofstream Hist, DATA, POS;

 public:
  MD(std::string DIRECTORY, size_t run_number);
  MD(std::string DIRECTORY, size_t run_number, bool QUENCH_F);
  ~MD();

  void Simulation(double DENSITY, double TEMPERATURE, int POWER, double A_CST);
  std::string get_dir();
  void reset_values();

 protected:
  void initialise(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                  std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz,
                  double TEMPERATURE);
  void mb_distribution(double TEMPERATURE, bool run_python_script);
  void verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry, std::vector<double> &rz,
                       std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz,
                       std::vector<double> &rrx, std::vector<double> &rry, std::vector<double> &rrz);
  void velocity_autocorrelation_function(std::vector<double> &Cvx, std::vector<double> &Cvy, std::vector<double> &Cvz);
  void radial_distribution_function(bool normalise = true);
  void mean_square_displacement(std::vector<double> &MSDx, std::vector<double> &MSDy, std::vector<double> &MSDz);
  void density_compression(int steps_quench, double TEMPERATURE);

  void open_files();
  void file_naming(int POWER, double A_cst);
  void write_to_files();
  void show_run(size_t step_size_show);
  void time_stamp(std::ofstream &, std::string variables);
  std::string convert_to_string(const double &x, const int &precision);
	std::string find_and_replace(std::string& source, const std::string& find, const std::string& replace);
};