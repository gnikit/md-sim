#include "MD.h"
#include <chrono>  // CPU run-time
#include <cstdint>
#include <ctime>    // std::chrono
#include <iomanip>  // setprecision
#include <iostream>
#include <iterator>
#include <sstream>        // stringstream
#include "FileLoading.h"  // FileLoading class

// Detects compiler and uses appropriate library
#if defined(__INTEL_COMPILER)
#include <mathimf.h>  // Intel Math library
#define COMPILER "INTEL"
#elif defined(__GNUC__)
#include <math.h>
#define COMPILER "G++"
#else
#include <math.h>
#define COMPILER "OTHER COMPILER"
#endif

// if changed, new vx,vy,vz files need to be generated
#define PARTICLES_PER_AXIS 10
// TODO: Increase the number of bins to find RDF intersects
#define NHIST 300                // Number of histogram bins
#pragma warning(disable : 4996)  //_CRT_SECURE_NO_WARNINGS

// Detects the OS and fetches the executable path that is passed
// in the FileLoadding.h class
#ifdef _WIN32
#define _WIN32 _WIN32
#include <windows.h>
std::string getExePath() {
  char result[MAX_PATH];
  return std::string(result, GetModuleFileName(NULL, result, MAX_PATH));
}
#else
#define _WIN32 0
#include <linux/limits.h>
#include <unistd.h>
std::string getExePath() {
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  return std::string(result, (count > 0) ? count : 0);
}
#endif

// TODO: Boltzmann Dist normalisation of the particles velocities in the
// beggining make it C++
// TODO: Calls to Python script are not multithread safe, calling a Py script
//       and reading from a file from multiple process could potentially lead to
//       undefined behaviour. Further testing required to prove this is not an
//       issue.

MD::MD(std::string DIRECTORY, size_t run_number) {
  _dir = DIRECTORY;
  _STEPS = run_number;

  Nx = Ny = Nz = PARTICLES_PER_AXIS;  // Number of particles per axis
  N = Nx * Ny * Nz;

  gr.resize(NHIST + 1, 0);  // gr with Index igr
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
  Cr.reserve(_STEPS);
  msd.reserve(_STEPS);
  u_en.reserve(_STEPS);
  k_en.reserve(_STEPS);
  pc.reserve(_STEPS);
  pk.reserve(_STEPS);
  temperature.reserve(_STEPS);
  PI = acos(-1.0);
}

MD::MD(std::string DIRECTORY, size_t run_number, bool QUENCH_F) {
  _dir = DIRECTORY;
  _STEPS = run_number;

  Nx = Ny = Nz = PARTICLES_PER_AXIS;  // Number of particles per axis
  N = Nx * Ny * Nz;

  gr.resize(NHIST + 1, 0);  // gr with Index igr
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
  Cr.reserve(_STEPS);
  msd.reserve(_STEPS);
  u_en.reserve(_STEPS);
  k_en.reserve(_STEPS);
  pc.reserve(_STEPS);
  pk.reserve(_STEPS);
  temperature.reserve(_STEPS);

  compression_flag = QUENCH_F;
  PI = acos(-1.0);
}
MD::~MD() {}

// Methods for MD Analysis
void MD::initialise(std::vector<double> &x, std::vector<double> &y,
                    std::vector<double> &z, std::vector<double> &vx,
                    std::vector<double> &vy, std::vector<double> &vz,
                    double TEMPERATURE) {
  /*
   *  Initialises the:
   *  + Position Arrays
   *  + Velocity Arrays (assign random velocities)
   *  + Conserves/ Scales momentum == 0
   *  + Temperature
   *  + Velocity Autocorrelaion Function
   *
   *  @param &x, &y, &z: X, Y, Z vector points
   *  @param &vx, &vy, &vz: Vx, Vy, Vz vector points
   *  @param TEMPERATURE: Thermostat target temperature
   */
  // Initialise position matrix and velocity matrix from Cubic Centred Lattice
  // TODO: Initialise from last position, FIX
  if (compression_flag == false) {
    size_t n = 0;
    size_t i, j, k;
    for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
        for (k = 0; k < Nz; k++) {
          x.push_back((i + 0.5) * scale);
          y.push_back((j + 0.5) * scale);
          z.push_back((k + 0.5) * scale);

          rrx.push_back((i + 0.5) * scale);
          rry.push_back((j + 0.5) * scale);
          rrz.push_back((k + 0.5) * scale);

          ++n;
        }
      }
    }
    // Generates Maxwell-Boltzmann dist from Python script
    mb_distribution(TEMPERATURE, true);
  }

  if (compression_flag == true && Q_counter == 0) {
    FileLoading<double> load_data;
    get_dir();  // initialises top_exe_dir
    std::string file_name = top_exe_dir +
                            "/data/Positions_Velocities_particles_" +
                            std::to_string(N) + ".txt";
    std::cout << "Try and read file: " << file_name << std::endl;
    std::vector<std::vector<double>> vel = load_data.LoadTxt(file_name, 9, '#');
    x = vel[0];
    y = vel[1];
    z = vel[2];
    vx = vel[3];
    vy = vel[4];
    vz = vel[5];
    rrx = vel[0];
    rry = vel[1];
    rrz = vel[2];
    std::cout << "Read successful" << std::endl;
    // Start from a highly thermalised fluid state
    // Generation of velocities with T = 10
    // MBDistribution(10);
    // TODO: Not sure what follows is correct in if-loop
    // Temperature calculation for the first step with very high T
    // scale of x, y, z
  }

  // scale of x, y, z
  double mean_vx = 0;
  double mean_vy = 0;
  double mean_vz = 0;

  size_t i;
  // Momentum conservation
  for (i = 0; i < N; i++) {
    mean_vx += vx[i] / N;
    mean_vy += vy[i] / N;
    mean_vz += vz[i] / N;
  }

  size_t tempN = N;
#pragma parallel
#pragma loop count min(128)
  // Subtracting Av. velocities from each particle
  for (i = 0; i < tempN; i++) {
    vx[i] = vx[i] - mean_vx;
    vy[i] = vy[i] - mean_vy;
    vz[i] = vz[i] - mean_vz;
  }
  // Temperature calculation, statistically
  KE = 0;
  for (i = 0; i < N; i++) {
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  T = KE / (1.5 * N);
  scale_v = sqrt(TEMPERATURE / T);  // scalling factor

  // Velocity scaling
#pragma parallel
#pragma loop count min(128)
  for (i = 0; i < tempN; i++) {
    vx[i] *= scale_v;
    vy[i] *= scale_v;
    vz[i] *= scale_v;
  }
  // MSD initialasation, storing first positions of particles
  MSDx = x;
  MSDy = y;
  MSDz = z;

  // VAF initialasation, storing first velocities of particles
  Cvx = vx;
  Cvy = vy;
  Cvz = vz;
  double first_val = 0;
  for (i = 0; i < N; i++) {
    first_val += (Cvx[i] * Cvx[i] + Cvy[i] * Cvy[i] + Cvz[i] * Cvz[i]) / N;
  }
  first_val /= N;
  Cr.push_back(first_val);
}

std::string MD::get_dir() {
  /*
   *	Returns the absolute, top working direcory of the git repo.
   *  Also, when called, the full path of the executable, along with its name
   *	will be stored in the full_exe_dir string.
   */
  // TODO: this is obvious duplication but see what happens to memory of
  // full_exe_path
  //      if nothing happens remove str
  std::string str = getExePath();
  full_exe_dir = getExePath();
  if (_WIN32) {
    // takes care of the windows backslashes
    full_exe_dir = find_and_replace(str, "\\", "/");
  }
  size_t stride = full_exe_dir.rfind("/MD-simulation");
  top_exe_dir = full_exe_dir.substr(0, stride) + "/MD-simulation";
  std::string rel_path = full_exe_dir.substr(stride + 1);

  std::cout << "Full executable directory: " << full_exe_dir << std::endl;
  std::cout << "Absolute top level directory: " << top_exe_dir << std::endl;

  return top_exe_dir;
}

void MD::mb_distribution(double TEMPERATURE, bool run_python_script = false) {
  std::string t = convert_to_string(TEMPERATURE, 4);
  std::string particles = std::to_string(N);
  std::string dir_str = get_dir();
  dir_str += "/data";

  if (run_python_script) {
    // Could be stored as variables and passed into FileNaming
    // rather than repeating the process
    // store in _particles_to_str, _T_to_str
    // taking care of the idiots that use spaces in paths
    std::string command =
        "python \"" + dir_str + "/MBDistribution.py\" " + particles + " " + t;
    system(command.c_str());  // Creates files with MD velocities
  }

  std::string vel_id = "_particles_" + particles + "_T_" + t + ".txt";
  FileLoading<double> obj;
  vx = obj.LoadSingleCol(dir_str + "/vx" + vel_id);
  vy = obj.LoadSingleCol(dir_str + "/vy" + vel_id);
  vz = obj.LoadSingleCol(dir_str + "/vz" + vel_id);
  std::cout << "Files loaded successfuly." << std::endl;
  // TODO: define in heap and delete FileLoading obj
}

void MD::verlet_algorithm(std::vector<double> &rx, std::vector<double> &ry,
                          std::vector<double> &rz, std::vector<double> &vx,
                          std::vector<double> &vy, std::vector<double> &vz,
                          std::vector<double> &rrx, std::vector<double> &rry,
                          std::vector<double> &rrz) {
  /*
   *  An iterative leap-frog Verlet Algorithm
   *
   *  @params <r> (x,y,z): position vector of particles.
   *  @params <v> <vx,vy,vz): velocity vector of particles.
   *  @params <rr> (xx,yy,zz): position vectors for the fluid in the
   *                           square displacement arrays.
   */
  size_t i;
  for (i = 0; i < N; i++) {
    vx[i] = vx[i] * scale_v + fx[i] * dt;
    vy[i] = vy[i] * scale_v + fy[i] * dt;
    vz[i] = vz[i] * scale_v + fz[i] * dt;
    rx[i] = rx[i] + vx[i] * dt;
    ry[i] = ry[i] + vy[i] * dt;
    rz[i] = rz[i] + vz[i] * dt;

    rrx[i] = rrx[i] + vx[i] * dt;
    rry[i] = rry[i] + vy[i] * dt;
    rrz[i] = rrz[i] + vz[i] * dt;

    // Kinetic Energy Calculation
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    // Boundary conditions Updated
    if (rx[i] > L) {
      rx[i] = rx[i] - L;
    }
    if (ry[i] > L) {
      ry[i] = ry[i] - L;
    }
    if (rz[i] > L) {
      rz[i] = rz[i] - L;
    }
    if (rx[i] < 0.0) {
      rx[i] = rx[i] + L;
    }
    if (ry[i] < 0.0) {
      ry[i] = ry[i] + L;
    }
    if (rz[i] < 0.0) {
      rz[i] = rz[i] + L;
    }
  }
}

void MD::velocity_autocorrelation_function(std::vector<double> &Cvx,
                                           std::vector<double> &Cvy,
                                           std::vector<double> &Cvz) {
  double temp = 0;  // resets every time step
  size_t i;
  for (i = 0; i < N; i++) {
    temp += (Cvx[i] * vx[i] + Cvy[i] * vy[i] + Cvz[i] * vz[i]);
  }
  temp /= N;
  Cr.push_back(temp);
}

void MD::radial_distribution_function(bool normalise) {
  /*normalise by default is TRUE*/
  // TODO: save raw RDF and normalisation factors & distance!
  double R = 0;
  double norm = 1;
  double cor_rho = _rho * (N - 1) / N;
  double dr = rg / NHIST;
  size_t i;
  Hist << "#particles (N): " << N << " steps: " << _STEPS << " rho: " << _rho << " bin (NHIST): " << NHIST << " cut_off (rg): "
       << rg << " dr: " << dr << std::endl;
  Hist << "#Unormalised" << '\t' << "Normalised" << std::endl;
  for (i = 1; i < NHIST; i++) {  // Changed initial loop value from 0 -> 1
    if (normalise) {
      R = rg * i / NHIST;
      // norm = (cor_rho * 2 * PI * R * R * N * _STEPS * dr);
      norm = cor_rho * (2.0 * PI * (pow((R + dr / 2.0), 3) - pow((R - dr / 2.0), 3)));
    }
    // gr[i] /= norm;  // not really needed
    Hist << gr[i] << std::endl;
  }
}

void MD::mean_square_displacement(std::vector<double> &MSDx,
                                  std::vector<double> &MSDy,
                                  std::vector<double> &MSDz) {
  double msd_temp = 0;
  for (size_t i = 0; i < N; ++i) {
    msd_temp += (pow((rrx[i] - MSDx[i]), 2) + pow((rry[i] - MSDy[i]), 2) +
                 pow((rrz[i] - MSDz[i]), 2));
  }
  msd_temp /= N;
  msd.push_back(msd_temp);
}

void MD::density_compression(int steps_quench, double TEMPERATURE) {
  // Increase _rho by 0.01
  _rho += 0.01;
  // Re-using this piece of code from MD::Simulation
  scale = pow((N / _rho), (1.0 / 3.0)) / PARTICLES_PER_AXIS;
  L = pow((N / _rho), 1.0 / 3.0);
  Vol = N / _rho;
  initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);
}

// MD Simulation
void MD::Simulation(double DENSITY, double TEMPERATURE, int POWER,
                    double A_CST) {
  /*
   *  Executes the fluid simulation. It includes all statistical methods
   *  defined in this class. It monitors the following quantities
   *  RDF, MSD, VAF, average particle energies & pressures.
   *
   *  @param DENSITY: the density rho of fluid.
   *  @param TEMPERATURE: the temperature of the thermostat.
   *  @param POWER: the power n that the pair potential will be using
   *                typical values are in the range of 6-18.
   *  @param A_CST: softening constant 'a' of the pair potential.
   *                When 'a' = 0, then fluid is pure MD, increasing
   *                'a' reuslts into softening of the pair potential.
   */
  // Initialise scalling variables
  std::cout << "***************************\n"
               "** MD Simulation started **\n"
               "***************************\n"
            << std::endl;
  std::cout << "Fluid parameters: rho: " << DENSITY << " T: " << TEMPERATURE
            << " n: " << POWER << " a: " << A_CST << std::endl;
  _T0 = TEMPERATURE;
  _rho = DENSITY;
  dt /= sqrt(_T0);
  // Box length scalling
  scale = pow((N / _rho), (1.0 / 3.0)) / PARTICLES_PER_AXIS;
  L = pow((N / _rho), 1.0 / 3.0);
  Vol = N / _rho;

  // cut_off redefinition
  // Large cut offs increase the runtime exponentially
  cut_off = 3.0;  // L / 2.;
  rg = cut_off;
  dr = rg / NHIST;

  // Filenaming should not be called
  file_naming(POWER, A_CST);
  open_files();
  time_stamp(DATA, "# step \t rho \t U \t K \t Pc \t Pk \t MSD \t VAF");

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);

  double xx, yy, zz;
  for (_STEP_INDEX = 0; _STEP_INDEX < _STEPS; _STEP_INDEX++) {
    // Forces loop
    // Resetting forces
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    U = 0;  // seting Potential U to 0
    PC = 0;

    size_t steps_quench = 10;  // steps between each quenching
    if (compression_flag == true && _STEP_INDEX != 0 &&
        _STEP_INDEX % steps_quench == 0) {
      ++Q_counter;
      density_compression(steps_quench, TEMPERATURE);
      std::cout << "Compressing fluid, run: " << Q_counter << " rho: " << _rho
                << std::endl;
    }

    size_t i, j;
    for (i = 0; i < N - 1; i++) {
      for (j = i + 1; j < N; j++) {
        x = rx[i] - rx[j];  // Separation distance
        y = ry[i] - ry[j];  // between particles i and j
        z = rz[i] - rz[j];  // in Cartesian

        xx = x;
        yy = y;
        zz = z;

        // Transposing elements with Periodic BC
        // xx, yy, zz are for the MSD calculation
        if (x > (0.5 * L)) {
          x = x - L;
          xx = xx - L;
        }
        if (y > (0.5 * L)) {
          y = y - L;
          yy = yy - L;
        }
        if (z > (0.5 * L)) {
          z = z - L;
          zz = zz - L;
        }
        if (x < (-0.5 * L)) {
          x = x + L;
          xx = xx + L;
        }
        if (y < (-0.5 * L)) {
          y = y + L;
          yy = yy + L;
        }
        if (z < (-0.5 * L)) {
          z = z + L;
          zz = zz + L;
        }

        r = sqrt((x * x) + (y * y) + (z * z));
        // TODO: enable q for BIP potential
        long double q = sqrt(r * r + A_CST * A_CST);

        // Force loop
        if (r < cut_off) {
          // TODO: implement functionally different potentials, currently
          //		using comment-uncomment to implement
          // BIP potential of the form: phi = 1/[(r**2 + a**2)**(n/2)]
          // TODO: BIP force
          long double ff =
              (POWER)*r * pow(q, ((-POWER - 2.0)));  // Force for particles

          // TODO: Gausian-force with sigma=1 and epsilon=1
          // long double ff = 2 * r * exp(-r * r);

          // Canceling the ij and ji pairs
          // Taking the lower triangular matrix
          fx[i] += x * ff / r;
          fx[j] -= x * ff / r;
          fy[i] += y * ff / r;
          fy[j] -= y * ff / r;
          fz[i] += z * ff / r;
          fz[j] -= z * ff / r;

          PC += r * ff;
          // TODO:Gaussian-Potential configurational Pressure
          // integral not evaluated

          // TODO: Add infinity and edge correction, do same for Pc

          // TODO: BIP potential
          U += pow(q, (-POWER));

          // TODO: Gaussian Potential GCM
          // U += exp(-r * r);

          // Radial Distribution
          igr = round(NHIST * r / rg);
          gr[igr] += 1;
          // rn = (igr - 0.5)*dr;
        }
      }
    }

    // Average Potential Energy per particle
    u_en.push_back(U / N);

    // Average Configurational Pressure Pc
    pc.push_back(PC / (3 * Vol));

    // Isothermal Calibration
    scale_v = sqrt(_T0 / T);  // using T & KE from prev timestep
    KE = 0;                   // resetting Kintetic Energy per iteration

    verlet_algorithm(rx, ry, rz, vx, vy, vz, rrx, rry, rrz);

    mean_square_displacement(MSDx, MSDy, MSDz);

    velocity_autocorrelation_function(Cvx, Cvy, Cvz);

    // Average Temperature
    T = KE / (1.5 * N);
    temperature.push_back(T);

    // Kinetic Pressure
    PK = _rho * T;
    pk.push_back(PK);

    // Average Kintetic Energy
    KE /= N;
    k_en.push_back(KE);

    // Density
    density.push_back(_rho);

    // ShowRun(500);  // shows every 500 steps
  }
  // Simulation Ends HERE

  write_to_files();
  // Saving Last Position
  time_stamp(POS, "# X\tY\tZ\tVx\tVy\tVz\tFx\tFy\tFz");
  for (size_t el = 0; el < rx.size(); el++) {
    POS << rx[el] << '\t' << ry[el] << '\t' << rz[el] << '\t' << vx[el] << '\t'
        << vy[el] << '\t' << vz[el] << '\t' << fx[el] << '\t' << fy[el] << '\t'
        << fz[el] << std::endl;
  }

  radial_distribution_function(true);  // normalisation argument
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
               "** MD Simulation terminated **\n"
               "******************************\n"
            << std::endl;
  // Streams should not close, vectors should not be cleared if object is to be
  // reused
  reset_values();  // no need to call if object is not reused
}

// File Handling
void MD::file_naming(int POWER, double A_cst) {
  /*
   * Generates a unique filename for the simulation results to be stored.
   * Three filenames are created for the RDF, last particles' states (position
   * ,velocity, acceleration) and statistical data.
   */
  std::stringstream A_stream, rho_stream, T_stream;

  T_stream << std::fixed << std::setprecision(4) << _T0;     // 4 decimal
  A_stream << std::fixed << std::setprecision(5) << A_cst;   // 5 decimals
  rho_stream << std::fixed << std::setprecision(4) << _rho;  // 4 decimal

  _step_to_str = "_step_" + std::to_string(_STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();
  _n_to_str = "_n_" + std::to_string(POWER);
  _A_to_str = "_A_" + A_stream.str();

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str + _T_to_str +
             _n_to_str + _A_to_str;

  // Explicit defitions
  _FILE_EXT = ".txt";
  data = "Data";
  pos = "Positions_Velocities";
  HIST = "RDF";

  // Path addition
  data = _dir + data + _FILE_ID + _FILE_EXT;
  pos = _dir + pos + _FILE_ID + _FILE_EXT;
  HIST = _dir + HIST + _FILE_ID + _FILE_EXT;
}

void MD::open_files() {
  /*
   * Open/Create if file does not exist
   * Overwrite existing data
   */
  Hist.open(HIST, std::ios::out | std::ios::trunc);
  DATA.open(data, std::ios::out | std::ios::trunc);
  POS.open(pos, std::ios::out | std::ios::trunc);
}

void MD::write_to_files() {
  /*
   * Writes values of parameters to file
   */
  for (size_t i = 0; i < _STEPS; i++) {
    DATA << (i + 1) << '\t' << density[i] << '\t' << temperature[i] << '\t'
         << u_en[i] << '\t' << k_en[i] << '\t' << pc[i] << '\t' << pk[i] << '\t'
         << msd[i] << '\t' << Cr[i] << std::endl;
  }
}

void MD::show_run(size_t step_size_show) {
  /*
   * Displays the system parameters every step_size_show of steps
   * Input the increment step
   */
  if (_STEP_INDEX == 0) {
    std::cout << "step:\tT:\tKE:\tU:\tU+K:\tPC:\tPK:\t(PK+PC):" << std::endl;
  }

  if (_STEP_INDEX % step_size_show == 0 || _STEP_INDEX == 1) {
    std::cout.precision(5);
    std::cout << _STEP_INDEX << "\t" << T << "\t" << KE << "\t" << U << "\t"
              << (U + KE) << "\t" << PC << "\t" << PK << "\t" << (PK + PC)
              << std::endl;
  }
}

void MD::reset_values() {
  /*
   * Closes open file streams and resets sizes and values to 0
   * For multiple simulations
   */
  // Close streams
  Hist.close();
  DATA.close();
  POS.close();

  // Clear values, size, but reserve capacity
  rx.clear();
  ry.clear();
  rz.clear();
  rrx.clear();
  rry.clear();
  rrz.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  density.clear();
  temperature.clear();
  u_en.clear();
  k_en.clear();
  pc.clear();
  pk.clear();
  msd.clear();
  Cr.clear();
  gr.resize(NHIST + 1, 0);  // gr with Index igr
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
}

void MD::time_stamp(std::ofstream &stream, std::string variables) {
  /*
   * Dates the file and allows the input of a header
   * Input a file stream to write and string of characters to display as headers
   */
  std::chrono::time_point<std::chrono::system_clock> instance;
  instance = std::chrono::system_clock::now();
  std::time_t date_time = std::chrono::system_clock::to_time_t(instance);
  stream << "# Created on: " << std::ctime(&date_time);
  stream << variables << std::endl;
}

std::string MD::convert_to_string(const double &x, const int &precision) {
  std::ostringstream ss;
  ss.str(std::string());  // don't forget to empty the stream
  ss << std::fixed << std::setprecision(precision) << x;

  return ss.str();
}

std::string MD::find_and_replace(std::string &source, const std::string &find,
                                 const std::string &replace) {
  for (std::string::size_type i = 0;
       (i = source.find(find, i)) != std::string::npos;) {
    source.replace(i, find.length(), replace);
    i += replace.length();
  }
  return source;
}
