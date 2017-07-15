//////////////////////////////////////////////////////////////////////
// Ioannis Nikiteas 13/7/2017									    //
//																    //
// BSc Dissertation:											    //
//  Investigating the transition from Molecular Dynamics to		    //
//	Smoothed Particle Hydrodynamics					    		    //
//																    //
//	University: Royal Holloway University of London		            //
//																    //
//	A program meant to simulate a MD fluid with an only	            //
//	repulsive pair-potential. Increasing the parameter A		 	//
//	creates a coarse-graining effect for the system allowing it		//
//	to transition to SPH											//
//																	//
//																	//
//	Use Intel C++ compiler with with MPI libraries for				//
//	optimal results.												//	   
//																	//	   
//////////////////////////////////////////////////////////////////////
#pragma once
#include <iostream>
#include <cmath>
#include <iomanip>  // setprecision
#include <vector>
#include <chrono>   // CPU run-time
#include <ctime>
#include <fstream>  // file writing
#include <iterator>
#include <assert.h>
#include <sstream>




class MD
{

protected:
	typedef std::vector<double> vec1d;	// short-hand notation, use with care

	vec1d rx;	vec1d ry;	vec1d rz;	// Position Arrays
	vec1d vx;	vec1d vy;	vec1d vz;	// Velocity Arrays
	vec1d fx;	vec1d fy;	vec1d fz;
	vec1d Cvx;	vec1d Cvy;	vec1d Cvz;
	vec1d rrx;	vec1d rry;	vec1d rrz;	// used in MSD
	vec1d MSDx;	vec1d MSDy;	vec1d MSDz;
	vec1d Cr;	vec1d msd;								// correlation vector with time index

	size_t Nx;	size_t Ny;	size_t Nz;
	size_t N; 	size_t N_step; size_t N_max;
	int power; double A;
	double dt = 0.005;		// time step dt = 0.005/sqrt(T0)
	double x, y, z;			// distance between particle i and j
	double r;				// distance in polar
	double rho;				// density
	double scale;
	double KE = 0.0;		// Kinetic Energy
	double T;				// Temperature
	double T0 = 1.4;		// Target Temperature. Desired T for the system to operate
	double L;				// Length of the box after scaling
							//std::cout << "Box side length: " << L << std::endl;
	double Vol;
	double cut_off = 2.5;	// Cut off radius for the analysis
	double U = 0;			// Potential Energy
	double PC = 0;			// Config Pressure
	double PK;
	double scale_v;

	// HISTOGRAM VARIABLES
	int Nhist, igr;		// Index of Hist
	double rg, den, rn;
	double dr;
	vec1d gr;

private:
	const long double pi = acos(-1.0);
	std::string run; std::string A_par; std::string separator;
	std::string path;
	std::string file_type;
	std::string HIST; std::string _VAF; std::string _MSD; std::string data;
	std::string pos;
	std::string Ldrx, Ldry, Ldrz, Ldvx, Ldvy, Ldvz;
	std::string _dir, _density, _step;
	std::ofstream Hist, VAF, MSD, DATA, POS;

public:

	MD(std::string DIRECTORY, double DENSITY, size_t run_number);
	~MD();

	void Simulation(int POWER, double A_cst);

protected:
	void Initialise(vec1d &x, vec1d &y, vec1d &z,
					vec1d &vx, vec1d &vy, vec1d &vz);
	void VerletAlgorithm(vec1d &rx, vec1d &ry, vec1d &rz,
						 vec1d &vx, vec1d &vy, vec1d &vz,
						 vec1d &rrx, vec1d &rry, vec1d &rrz);
	void VelocityAutocorrelationFunction(vec1d &Cvx, vec1d &Cvy, vec1d &Cvz);
	void RadialDistributionFunction();
	void MeanSquareDisplacement(vec1d &MSDx, vec1d &MSDy, vec1d &MSDz);

	void OpenFiles();
	void CreateFiles(int POWER, double A_cst);
	void WriteToFiles();
	void ShowRun(size_t step_size_show);
	void ResetValues();
	void time(std::ofstream&, std::string variables);
};
