//////////////////////////////////////////////////////////////////////
// Ioannis Nikiteas 13/7/2017									                      //
//																                                  //
// BSc Dissertation:										                      	    //
//  Investigating the transition from Molecular Dynamics to		      //
//	Smoothed Particle Hydrodynamics					    		                //
//																                                  //
//	University: Royal Holloway University of London		              //
//																                                  // 
//	A program meant to simulate a MD fluid with an only	            //
//	repulsive pair-potential. Increasing the parameter A		 	      //
//	creates a coarse-graining effect for the system allowing it		  //
//	to transition to SPH											                      //
//																	                                //
//																	                                //
//																	                                //
//////////////////////////////////////////////////////////////////////
#pragma once
#include <iostream>
#include <mathimf.h>
#include <iomanip>  // setprecision
#include <vector>
#include <chrono>   // CPU run-time
#include <ctime>
#include <fstream>  // file writing
#include <iterator>
#include <assert.h>
#include <sstream>
#include <cstdint>

class MD {
protected:
  typedef std::vector<double> vec1d;	// short-hand notation, use with care

  vec1d rx, ry, rz;	// Position Arrays
  vec1d vx, vy, vz;	// Velocity Arrays
  vec1d fx, fy, fz;
  vec1d Cvx, Cvy, Cvz;
  vec1d rrx, rry, rrz;	// used in MSD
  vec1d MSDx, MSDy, MSDz;
  vec1d Cr, msd;								// correlation vector with time index

  size_t Nx, Ny, Nz;   // Particles in x, y, z
  size_t N, _STEP_INDEX, _STEPS;  // Total particles, step index, maximum steps
  double _T0;		// Target Temperature. Desired T for the system to operate
  double dt = 0.005;	// time step dt = 0.005/sqrt(T0)
  double x, y, z;			// distance between particle i and j
  double r;				    // distance in polar
  double _rho;				// density
  double scale;
  double KE = 0.0;		// Kinetic Energy
  double T;				    // Temperature

  double L;				    // Length of the box after scaling
  double Vol;
  double cut_off = 2.5;	// Cut off radius for the analysis
  double U = 0;			    // Potential Energy
  double PC = 0;			  // Config Pressure
  double PK;
  double scale_v;

  // HISTOGRAM VARIABLES
  int igr;		// Index of Hist
  double rg;
  double dr;
  vec1d gr;

private:
  const long double PI = acos(-1.0);
  std::string _FILE_EXT;
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str, _n_to_str, _A_to_str;
  std::string HIST, _VAF, _MSD, data, pos; // contain full names with args
  std::string _dir, _FILE_ID;
  std::ofstream Hist, VAF, MSD, DATA, POS;

public:

  MD(std::string DIRECTORY, size_t run_number);
  ~MD();

  void Simulation(double DENSITY, double TEMPERATURE, int POWER, double A_CST);
  std::string getDir();

protected:
  void Initialise(vec1d &x, vec1d &y, vec1d &z,
                  vec1d &vx, vec1d &vy, vec1d &vz,
                  double TEMPERATURE);
  void MBDistribution(double TEMPERATURE);
  void VerletAlgorithm(vec1d &rx, vec1d &ry, vec1d &rz,
                       vec1d &vx, vec1d &vy, vec1d &vz,
                       vec1d &rrx, vec1d &rry, vec1d &rrz);
  void VelocityAutocorrelationFunction(vec1d &Cvx, vec1d &Cvy, vec1d &Cvz);
  void RadialDistributionFunction();
  void MeanSquareDisplacement(vec1d &MSDx, vec1d &MSDy, vec1d &MSDz);

  void OpenFiles();
  void FileNaming(int POWER, double A_cst);
  void WriteToFiles();
  void ShowRun(size_t step_size_show);
  void ResetValues();
  void TimeStamp(std::ofstream&, std::string variables);
  std::vector<double> ReadFromFile(const std::string &file_name);
  std::string ConvertToString(const double & x, const int & precision);
};