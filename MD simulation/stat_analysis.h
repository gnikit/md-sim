#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <fstream>	
#include <assert.h>
#include <iterator>
#include <sstream>
#include <iomanip>  // setprecision
#include <exception>

class Stat_Analysis {

protected:
  typedef std::vector<double> vec1d;
  std::ifstream _data_reader;
  size_t _N, _STEPS;  // Number of particles
  double _RHO, _T0;


private:

  vec1d T_vec, K_vec, U_vec, E_vec, Pc_vec, Pk_vec, P_vec;
  vec1d _A_list;
  std::string _path;
  std::string _FILE_EXT;
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str, _n_to_str, _A_to_str;
  std::string _DIR, _FILE_ID;

  long double _temp0 = 0;
  long double _temp1 = 0;
  long double _temp2 = 0;
  long double _temp3 = 0;
  long double _temp4 = 0;
  long double _temp5 = 0;
  long double _temp6 = 0;

public:
  Stat_Analysis(std::string PATH, vec1d A_LIST, size_t STEPS, size_t T, size_t PARTICLES, double DENSITY);
  ~Stat_Analysis();

  void ReadFromFile(vec1d &x, const std::string &file_name);
  void ReadFromFile(vec1d &T, vec1d &K, vec1d &U,
                    vec1d &E, vec1d &Pc, vec1d &Pk, vec1d &P,
                    const std::string &file_name);
  void Mean(vec1d &T, vec1d &K, vec1d &U,
            vec1d &E, vec1d &Pc, vec1d &Pk, vec1d &P);
  void StaticDataProcessing(size_t n);
  // TODO: make a RunMe method and add arguments to the existing methods for increased reusability

};