#pragma once
#include <chrono>
#include <fstream>
#include <iomanip>  // setprecision
#include <sstream>  // stringstream
#include <string>
#include <vector>

class Stat_File {
 public:
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str;
  std::string _n_to_str, _A_to_str;
  std::string rdf, data, pos;
  std::string _dir, _FILE_ID;
  std::ofstream RDF, DATA, POS;

  // Initialise them before doing any IO operations
  size_t STEPS, N;

  Stat_File();
  Stat_File(size_t &STEPS, size_t &N, std::string DIR);
  void open_files(std::string &data, std::string &rdf, std::string &pos);

  void write_data_file(std::vector<double> &density,
                       std::vector<double> &temperature,
                       std::vector<double> &u_en, std::vector<double> &k_en,
                       std::vector<double> &pc, std::vector<double> &pk,
                       std::vector<double> &msd, std::vector<double> &Cr,
                       std::vector<double> &sfx, std::vector<double> &sfy,
                       std::vector<double> &sfz);
  void time_stamp(std::ofstream &, std::string variables);
  std::string file_naming(std::string prefix, double &DENSITY,
                          double &TEMPERATURE, double &POWER, double &A_cst);
  std::string convert_to_string(const double &x, const int &precision);
};