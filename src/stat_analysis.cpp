#include "stat_analysis.h"

Stat_Analysis::Stat_Analysis(std::string PATH, size_t STEPS, size_t PARTICLES,
                             double DENSITY, double T, vec1d A_LIST) {
  _DIR = PATH;
  _A_list = A_LIST;
  _STEPS = STEPS;
  _T0 = T;
  _N = PARTICLES;
  _RHO = DENSITY;
}

Stat_Analysis::~Stat_Analysis() {}

void Stat_Analysis::ReadFromFile(vec1d &x, const std::string &file_name) {
  //TODO: this could be inherited from MD.cpp
  /*
	Reads from a stream that already exists for a file that is already placed in the
	directory and appends the data into a 1D vector.
	*/
  std::ifstream read_file(file_name);
  assert(read_file.is_open());

  std::copy(std::istream_iterator<long double>(read_file),
            std::istream_iterator<long double>(), std::back_inserter(x));

  read_file.close();
}

void Stat_Analysis::ReadFromFile(vec1d &T, vec1d &K, vec1d &U,
                                 vec1d &E, vec1d &Pc, vec1d &Pk, vec1d &P,
                                 const std::string &file_name) {
  /*
	  Reads from a stream that already exists for a file that is already placed in the
	  directory and appends the data into a 1D vector.
	*/
  std::ifstream read_file(file_name);
  try {
    if (read_file.is_open() == false) {
      throw "File Not Found";
    }
    // TODO: Fix exception throwing, currently exception is not caught properly.
    //assert(read_file.is_open()); // Evaluates false if file does not exist
    long double a, b, c, d, e, f, g;
    std::string line;
    while (!read_file.eof()) {  // Stops when End Of File is reached
      std::getline(read_file, line);

      //TODO: the IF and ELSE statements should be reversed and then ELSE should be removed
      if (line.length() == 0 || line[0] == '#') {
        //std::cout << "IGNORE LINE" << std::endl;
      } else {
        std::stringstream ss(line);
        ss >> a;
        ss >> b;
        ss >> c;
        ss >> d;
        ss >> e;

        T.push_back(a);
        K.push_back(b);
        U.push_back(c);
        E.push_back(d);
        Pc.push_back(e);
      }
    }
  } catch (std::string e) {
    std::cout << e << std::endl;
  }
  read_file.close();
}

void Stat_Analysis::Mean(vec1d &T, vec1d &K, vec1d &U,
                         vec1d &E, vec1d &Pc, vec1d &Pk, vec1d &P) {
  size_t n = 0;

  for (size_t i = 0; i < T.size(); i++) {
    _temp0 += T[i];
    _temp1 += K[i];
    _temp2 += U[i];
    _temp3 += E[i];
    _temp4 += Pc[i];
    //_temp5 += Zz[i];
    //_temp6 += w[i];
    n++;
  }
  _temp0 /= n;  // this is the mean
  _temp1 /= n;
  _temp2 /= n;
  _temp3 /= n;
  _temp4 /= n;
  //_temp5 /= n;
  //_temp6 /= n;
}

void Stat_Analysis::StaticDataProcessing(size_t POWER)
/*
Used to average the energies and pressures for all
values of A for a specific number n of the potential
strength
*/
{
  vec1d temp;

  std::stringstream A_stream;  // Fixing double to 2 decimals
  std::stringstream rho_stream;
  std::stringstream T_stream;
  // TODO: setprecission function input here

  T_stream << std::fixed << std::setprecision(2) << _T0;     // 1 decimal
  rho_stream << std::fixed << std::setprecision(2) << _RHO;  // 1 decimal

  // For consistency string streams should be placed here
  _step_to_str = "_step_" + std::to_string(_STEPS);        // constr
  _particles_to_str = "_particles_" + std::to_string(_N);  // constr
  _rho_to_str = "_rho_" + rho_stream.str();                // constr
  _T_to_str = "_T_" + T_stream.str();                      // constr
  _n_to_str = "_n_" + std::to_string(POWER);
  _A_to_str = "_A_" + A_stream.str();

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str +
             _T_to_str + _n_to_str;

  std::ofstream data;
  std::string _FILE_NAME, sep, power, a, TXT;
  sep = "_";
  TXT = ".txt";
  power = std::to_string(POWER);
  std::string name = _DIR + "AVGdata" + _FILE_ID + TXT;
  data.open(name, std::ios::out | std::ios::trunc);

  data << "# A:\tT\tK:\tU:\tETot:\tPc:" << std::endl;

  for (size_t i = 0; i < _A_list.size(); i++) {
    _FILE_NAME.clear();
    T_vec.clear();
    K_vec.clear();
    U_vec.clear();
    E_vec.clear();
    Pc_vec.clear();
    Pk_vec.clear();
    P_vec.clear();
    _FILE_NAME = "Data";

    // Converts A to 2 decimals
    std::stringstream stream;
    stream << std::fixed << std::setprecision(4) << _A_list[i];

    _A_to_str = "_A_" + stream.str();
    _FILE_NAME = _DIR + _FILE_NAME + _FILE_ID + _A_to_str + TXT;
    ReadFromFile(T_vec, K_vec, U_vec, E_vec, Pc_vec, Pk_vec, P_vec, _FILE_NAME);
    Mean(T_vec, K_vec, U_vec, E_vec, Pc_vec, Pk_vec, P_vec);
    data.precision(5);
    // File writting
    data << _A_list[i] << '\t' << _temp0 << '\t' << _temp1 << '\t' << _temp2 << '\t' << _temp3 << '\t' << _temp4 << std::endl;
    // consider excluding usless vectors like PK, PTOT
  }
}