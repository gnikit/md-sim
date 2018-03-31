#include "MD.h"
#define PARTICLES_PER_AXIS 10
#define NHIST 300
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
// TODO: Boltzmann Dist normalisation of the particles velocities in the beggining


MD::MD(std::string DIRECTORY, size_t run_number) {
  _dir = DIRECTORY;
  _STEPS = run_number;

  Nx = Ny = Nz = PARTICLES_PER_AXIS; // Number of particles per axis
  N = Nx * Ny * Nz;

  gr.resize(NHIST + 1, 0); // gr with Index igr
  fx.resize(N, 0); fy.resize(N, 0); fz.resize(N, 0);
  //TODO: reserve r's and v's, instead of push_back
}
MD::~MD() {}

// Methods for MD Analysis
void MD::Initialise(vec1d &x, vec1d &y, vec1d &z,
                    vec1d &vx, vec1d &vy, vec1d &vz, double TEMPERATURE) {
  /*
  * Initialises the:
  * + Position Arrays
  * + Velocity Arrays (assign random velocities)
  * + Conserves/ Scales momentum == 0
  * + Temperature
  * + Velocity Autocorrelaion Function
  */
  // Initialise position matrix and velocity matrix
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

        // TODO: Velocitied need to follow a Maxwell-Boltzmann Dist
        ++n;
      }
    }
  }
  // Reading Maxwell Boltzmann velocity Dist from files
  // TODO: Python script buggy with argument passing
  // directory defined wrt the dir wheere .o will execute
  // Windowd: PWD=MD-simulation
  // TODO: Add precompiler handles for Linux/Windows
  //       if Win -> load from /data/vx.txt etc.
  //       else ->  load from ../data/vx.txt etc. (makefile)
  vx = ReadFromFile("../data/vx.txt");
  vy = ReadFromFile("../data/vy.txt");
  vz = ReadFromFile("../data/vz.txt");
  //MBDistribution(TEMPERATURE);
  // scale of x, y, z
  double mean_vx = 0;
  double mean_vy = 0;
  double mean_vz = 0;

  // Momentum conservation array
  for (i = 0; i < N; i++) {
    mean_vx += vx[i] / N; // Calculating Average velocity for each dimension
    mean_vy += vy[i] / N;
    mean_vz += vz[i] / N;
  }

  size_t tempN = N;    // Recommended opt by Intel

  for (i = 0; i < tempN; i++) {
    vx[i] = vx[i] - mean_vx; // Subtracting Av. velocities from each particle
    vy[i] = vy[i] - mean_vy;
    vz[i] = vz[i] - mean_vz;
  }
  // T Calc
  KE = 0;
  for (i = 0; i < N; i++) {
    KE += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  T = KE / (1.5 * N);
  scale_v = sqrt(TEMPERATURE / T); // scalling factor

  // Velocity scaling
  for (i = 0; i < tempN; i++) {
    vx[i] *= scale_v;
    vy[i] *= scale_v;
    vz[i] *= scale_v;
  }
  // MSD initialasation
  MSDx = x;
  MSDy = y;
  MSDz = z;

  // VAF initialasation
  Cvx = vx;
  Cvy = vy;
  Cvz = vz;
  double first_val = 0;
  for (i = 0; i < N; i++) {
    first_val += (Cvx[i] * Cvx[i] + Cvy[i] * Cvy[i] + Cvz[i] * Cvz[i]) / N;
  }
  first_val /= N;
  //Cr.push_back(first_val);	// HACK: Enable for debugging
  VAF << first_val << std::endl;
}

void MD::MBDistribution(double TEMPERATURE) {
  std::string t = ConvertToString(TEMPERATURE, 2);
  std::string particles = std::to_string(N); // defined in constructor
  // Could be stored as variables and passed into FileNaming
  // rather than repeating the process
  // store in _particles_to_str, _T_to_str
  // TODO: Python script is buggy and sometimes gives error when passing arguments
  //std::string command = "python C:/Users/gn/source/repos/MD-simulation/MBDistribution.py " + particles + " " + t;
  //system(command.c_str());  // Creates files with MD velocities

  std::string vel_id = "particles_" + particles + "_T_" + t;

  vx = ReadFromFile(vel_id + "vx.txt");
  vy = ReadFromFile(vel_id + "vy.txt");
  vz = ReadFromFile(vel_id + "vz.txt");
}


void MD::VerletAlgorithm(vec1d &rx, vec1d &ry, vec1d &rz,
                         vec1d &vx, vec1d &vy, vec1d &vz,
                         vec1d &rrx, vec1d &rry, vec1d &rrz) {
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

void MD::VelocityAutocorrelationFunction(vec1d &Cvx,
                                         vec1d &Cvy,
                                         vec1d &Cvz) {
  double temp_ = 0; // resets every time step
  size_t i;
  for (i = 0; i < N; i++) {
    temp_ += (Cvx[i] * vx[i] + Cvy[i] * vy[i] + Cvz[i] * vz[i]);
  }
  temp_ /= N;
  //Cr.push_back(temp_);	// HACK: Enable for debugging
  VAF << temp_ << std::endl; // writes to file
}

void MD::RadialDistributionFunction() {
  double R = 0;
  double norm;
  double cor_rho = _rho * (N - 1) / N;
  size_t i;
  for (i = 1; i < NHIST; i++) {  // Changed initial loop value from 0 -> 1
    R = rg * i / NHIST;
    norm = (cor_rho * 2 * PI * R * R * N * _STEPS * dr);
    gr[i] /= norm;	// not really needed
    Hist << gr[i] << std::endl;
  }
}

void MD::MeanSquareDisplacement(vec1d &MSDx,
                                vec1d &MSDy,
                                vec1d &MSDz) {
  double msd_temp = 0;
  for (size_t i = 0; i < N; ++i) {
    msd_temp += (pow((rrx[i] - MSDx[i]), 2) + pow((rry[i] - MSDy[i]), 2) +
                 pow((rrz[i] - MSDz[i]), 2));
  }
  msd_temp /= N;
  //msd.push_back(msd_temp);	// HACK: Enable for debugging
  MSD << msd_temp << std::endl;
}

// MD Simulation
void MD::Simulation(double DENSITY, double TEMPERATURE, int POWER, double A_CST) {
  // Initialise scalling variables
  // If Simulation(...) is not run, _T0, _rho need to be initialised elsewhere
  _T0 = TEMPERATURE;
  _rho = DENSITY;
  dt /= sqrt(_T0);
  // Box length scalling
  scale = pow((N / _rho), (1.0 / 3.0)) / PARTICLES_PER_AXIS;
  L = pow((N / _rho), 1.0 / 3.0);
  Vol = N / _rho;

  //TASK: cut_off redefinition
  cut_off = 3.0;//L / 2.;
  rg = cut_off;
  dr = rg / NHIST;

  FileNaming(POWER, A_CST);
  OpenFiles();
  TimeStamp(DATA, "# T\tK\tU\tEtot\tPc\tPk\tPtot");
  std::chrono::steady_clock::time_point begin =
    std::chrono::steady_clock::now();
  Initialise(rx, ry, rz, vx, vy, vz, TEMPERATURE);

  double xx, yy, zz;
  for (_STEP_INDEX = 0; _STEP_INDEX < _STEPS; _STEP_INDEX++) {
    // Forces loop
    // Resetting forces
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    U = 0; // seting Potential U to 0
    PC = 0;
    size_t i, j;
    for (i = 0; i < N - 1; i++) {
      for (j = i + 1; j < N; j++) {
        x = rx[i] - rx[j]; // Separation distance
        y = ry[i] - ry[j]; // between particles i and j
        z = rz[i] - rz[j]; // in Cartesian

        xx = rrx[i] - rrx[j];
        yy = rry[i] - rry[j];
        zz = rrz[i] - rrz[j];

        // Transposing elements with Periodic BC
        if (x > (0.5 * L)) {
          x = x - L;
        }
        if (y > (0.5 * L)) {
          y = y - L;
        }
        if (z > (0.5 * L)) {
          z = z - L;
        }
        if (x < (-0.5 * L)) {
          x = x + L;
        }
        if (y < (-0.5 * L)) {
          y = y + L;
        }
        if (z < (-0.5 * L)) {
          z = z + L;
        }

        /////////// for MSD  /////////////////
        if (xx > (0.5 * L)) {
          xx = xx - L;
        }
        if (yy > (0.5 * L)) {
          yy = yy - L;
        }
        if (zz > (0.5 * L)) {
          zz = zz - L;
        }
        if (xx < (-0.5 * L)) {
          xx = xx + L;
        }
        if (yy < (-0.5 * L)) {
          yy = yy + L;
        }
        if (zz < (-0.5 * L)) {
          zz = zz + L;
        }

        r = sqrt((x * x) + (y * y) + (z * z));
        long double q = sqrt(r * r + A_CST * A_CST);

        // Force loop
        if (r < cut_off) {
          long double ff =
            (POWER)*r *	pow(q, ((-POWER - 2.0))); // Force for particles

          fx[i] += x * ff / r;
          fx[j] -= x * ff / r; // Canceling the ij and ji pairs
          fy[i] += y * ff / r; // Taking the lower triangular matrix
          fy[j] -= y * ff / r;
          fz[i] += z * ff / r;
          fz[j] -= z * ff / r;

          PC += r * ff;
          // TODO: Add infinity and edge correction, do same for Pc
          U += pow(q, (-POWER)); // Potential Calculation

          // Radial Distribution
          igr = round(NHIST * r / rg);
          gr[igr] += 1;
          //rn = (igr - 0.5)*dr;
        }
      }
    }

    U /= N; // Average Potential Energy per particle
    PC = PC / (3 * Vol);

    // Isothermal Calibration
    scale_v = sqrt(_T0 / T); // using T & KE from prev timestep
    KE = 0; // set 0 for each step

    // Verlet Algorithm
    VerletAlgorithm(rx, ry, rz, vx, vy, vz, rrx, rry, rrz);
    // MSD
    MeanSquareDisplacement(MSDx, MSDy, MSDz);
    // VAF
    VelocityAutocorrelationFunction(Cvx, Cvy, Cvz);

    T = KE / (1.5 * N); // Average T
    PK = _rho * T;       // Kinetic part of pressure
    KE /= N;

    WriteToFiles();
    //ShowRun(500);  // shows every 500 steps
  }
  // Simulation Ends HERE
  // Saving Last Position
  TimeStamp(POS, "# X\tY\tZ\tVx\tVy\tVz");
  for (size_t el = 0; el < rx.size(); el++) {
    POS << rx[el] << '\t' << ry[el] << '\t'
      << rz[el] << '\t' << vx[el] << '\t'
      << vy[el] << '\t' << vz[el] << std::endl;
  }

  RadialDistributionFunction();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
    << "CPU run time = "
    << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() /
    60
    << " min "
    << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() %
    60
    << "s" << std::endl;
  ResetValues(); // Check if everything is reset, could have missed something
}

std::string MD::getDir() {
  /*
  * Returns the directory in a string
  */
  return _dir;
}

// File Handling
void MD::FileNaming(int POWER, double A_cst) {
  /*
  * Generates file names for the different I/O operations
  */
  std::stringstream A_stream, rho_stream, T_stream;

  // TODO: setprecission function input here

  T_stream << std::fixed << std::setprecision(4) << _T0;    // 4 decimal
  A_stream << std::fixed << std::setprecision(5) << A_cst;  // 5 decimals
  rho_stream << std::fixed << std::setprecision(4) << _rho;	// 4 decimal

  _step_to_str = "_step_" + std::to_string(_STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();
  _n_to_str = "_n_" + std::to_string(POWER);
  _A_to_str = "_A_" + A_stream.str();

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str +
    _T_to_str + _n_to_str + _A_to_str;

  // Explicit defitions 
  _FILE_EXT = ".txt";
  data = "Data";
  pos = "Positions_Velocities";
  HIST = "Hist";
  _VAF = "VAF";
  _MSD = "MSD";

  // Path addition
  data = _dir + data + _FILE_ID + _FILE_EXT;
  pos = _dir + pos + _FILE_ID + _FILE_EXT;
  HIST = _dir + HIST + _FILE_ID + _FILE_EXT;
  _VAF = _dir + _VAF + _FILE_ID + _FILE_EXT;
  _MSD = _dir + _MSD + _FILE_ID + _FILE_EXT;
}

void MD::OpenFiles() {
  /*
  * Open/Create if file does not exist
  * Overwrite existing data
  */
  // opens for files output and deletes prev content
  Hist.open(HIST, std::ios::out | std::ios::trunc);
  VAF.open(_VAF, std::ios::out | std::ios::trunc);
  MSD.open(_MSD, std::ios::out | std::ios::trunc);
  DATA.open(data, std::ios::out | std::ios::trunc);
  POS.open(pos, std::ios::out | std::ios::trunc);
}

void MD::WriteToFiles() {
  /*
  * Writes values of parameters to file
  */
  DATA << T << '\t' << KE << '\t' << U << '\t'
    << (U + KE) << '\t' << PC << '\t' << PK
    << '\t' << (PC + PK) << std::endl;
}

void MD::ShowRun(size_t step_size_show) {
  /*
  *Displays the system parameters every step_size_show of steps
  *Input the increment step
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

void MD::ResetValues() {
  /*
  Closes open file streams and resets sizes and values to 0
  For multiple simulations
  */
  // Close streams at the end of run
  Hist.close();
  VAF.close();
  MSD.close();
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
  gr.resize(NHIST + 1, 0); // gr with Index igr
  fx.resize(N, 0);
  fy.resize(N, 0);
  fz.resize(N, 0);
}

void MD::TimeStamp(std::ofstream& stream, std::string variables) {
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

std::vector<double> MD::ReadFromFile(const std::string & file_name) {
  /*
  * Reads from a stream that already exists for a file that is already placed in the
  * directory and appends the data into a 1D vector.
  */
  std::vector<double> data;
  std::ifstream read_file(file_name);

  assert(read_file.is_open());

  std::copy(std::istream_iterator<long double>(read_file),
            std::istream_iterator<long double>(), std::back_inserter(data));

  read_file.close();
  return data;
}

std::string MD::ConvertToString(const double & x, const int & precision) {
  static std::ostringstream ss;
  ss.str(std::string()); // don't forget to empty the stream
  ss << std::fixed << std::setprecision(precision) << x;

  return ss.str();
}