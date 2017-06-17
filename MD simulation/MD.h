//////////////////////////////////////////////////////////////////////
// Ioannis Nikiteas 23/2/2017										//
//																	//
// BSc Dissertation:												//
//  Investigating the transition from Molecular Dynamics to			//
//	Smoothed Particle Hydrodynamics									//
//																	//
//	University: Royal Holloway University of London					//
//																	//
//	A program meant to simulate a MD fluid with an only				//
//	repulsive pair-potential. Increasing the parameter A			//
//	creates a coarse-graining effect for the system allowing it		//
//	to transition to SPH											//
//																	//
//																	//
//	Use Intel C++ compiler with O3 and GPU acceleration for			//
//	optimal results.												//	   
//																	//	   
//																	//
//																	//
//																	//
//																	//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <iostream>
#include <math.h>
#include <iomanip>  // setprecision
#include <vector>
#include <chrono>   // CPU run-time
#include <fstream>  // file writing
#include <iterator>
#include <assert.h>
#include <sstream>



class MD
{

protected:
	std::vector <double> rx;	std::vector <double> ry;	std::vector <double> rz;	// Position Arrays
	std::vector <double> vx;	std::vector <double> vy;	std::vector <double> vz;	// Velocity Arrays
	std::vector <double> fx;	std::vector <double> fy;	std::vector <double> fz;
	std::vector<double> Cvx;    std::vector<double> Cvy;    std::vector<double> Cvz;
	std::vector<double> rrx;	std::vector<double> rry;	std::vector<double> rrz;	// used in MSD
	std::vector<double> MSDx;	std::vector<double> MSDy;	std::vector<double> MSDz;
	std::vector<double> Cr;		std::vector<double> msd;								// correlation vector with time index


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
	size_t Nhist, igr;		// Index of Hist
	double rg, den, rn;
	double dr;
	std::vector<double> gr;

private:
	const long double pi = acos(-1.0);
	std::string run; std::string A_par; std::string separator;
	std::string path;
	std::string file_type;
	std::string KIN; std::string POT; std::string TOT;
	std::string PRESSUREC; std::string PRESSUREK; std::string PCKTOT; std::string TEMPERATURE;
	std::string HIST; std::string _VAF; std::string _MSD; std::string data;
	std::string Ldrx, Ldry, Ldrz, Ldvx, Ldvy, Ldvz;
	std::string _dir, _density, _step;
	std::ofstream Temperature, KinEn, PotEn, TotEn, PressureC, PressureK, PCK, Hist, VAF, MSD, DATA;
	std::ofstream Loadrx, Loadry, Loadrz, Loadvx, Loadvy, Loadvz;

public:

	MD(std::string DIRECTORY, double DENSITY, size_t run_number);
	~MD();

	void Simulation(int POWER, double A_cst);

protected:
	void Initialise(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
		std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz);
	void VerletAlgorithm(std::vector<double> &rx, std::vector<double> &ry, std::vector<double> &rz,
		std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz,
		std::vector<double> &rrx, std::vector<double> &rry, std::vector<double> &rrz);
	void VelocityAutocorrelationFunction(std::vector<double> &Cvx, std::vector<double> &Cvy, std::vector<double> &Cvz);
	void RadialDistributionFunction();
	void MeanSquareDisplacement(std::vector<double> &MSDx, std::vector<double> &MSDy, std::vector<double> &MSDz);

	void OpenFiles();
	void CreateFiles(int POWER, double A_cst);
	void WriteToFiles();
	void ShowRun(size_t step_size_show);
};
