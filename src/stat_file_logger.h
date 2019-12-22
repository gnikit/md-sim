#pragma once
#include <math.h>  // isnan
#include <chrono>
#include <fstream>
#include <iomanip>  // setprecision
#include <sstream>  // stringstream
#include <string>
#include <vector>

class stat_file {
 public:
  std::string _step_to_str, _particles_to_str, _rho_to_str, _T_to_str;
  std::string _n_to_str, _A_to_str;
  std::string _dir, _FILE_ID;
  std::vector<std::string> file_names;

  // Initialise them before doing any IO operations
  stat_file();

  static std::vector<std::ofstream> open_files(
      std::vector<std::string> const &file_names);

  static void write_data_file(
      std::ofstream &file_stream, std::string const &header,
      std::vector<std::vector<double>> const &all_output_vectors);

  static void time_stamp(std::ofstream &file_stream,
                         std::string const &variables);

  std::string file_naming(std::string const &prefix, size_t const &STEPS,
                          size_t const &N, double const &DENSITY,
                          double const &TEMPERATURE, double const &POWER,
                          double const &A_cst);

  static std::string convert_to_string(const double &x, const int &precision);

  void write_file(std::vector<std::vector<double>> &output_quantities,
                  std::ofstream &fstream, std::string const &header);
};