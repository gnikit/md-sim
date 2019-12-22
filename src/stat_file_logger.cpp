#include "stat_file_logger.h"
#include "FileIO.h"

// Passes the name of the variable to a stream
#define GET_VAR_NAME(stream, variable) (stream) << #variable

stat_file::stat_file() {}

std::vector<std::ofstream> stat_file::open_files(
    std::vector<std::string> const &file_names) {
  /**
   * Open/Create if file does not exist.
   * Overwrites existing data.
   */
  std::vector<std::ofstream> file_streams;
  for (size_t file = 0; file < file_names.size(); ++file) {
    std::ofstream temp;
    temp.open(file_names[file], std::ios::out | std::ios::trunc);
    file_streams.push_back(std::move(temp));
  }

  return file_streams;
}

void stat_file::write_data_file(
    std::ofstream &file_stream, std::string const &header,
    std::vector<std::vector<double>> const &all_output_vectors) {
  // Write the timestamp and header to the stream
  time_stamp(file_stream, header);

  // Find the largest size vector in all our vectors
  size_t rows = all_output_vectors[0].size();
  for (auto const &i : all_output_vectors)
    if (i.size() > rows) rows = i.size();
  
  for (size_t i = 0; i < rows; ++i) {
    std::string line = "";
    for (size_t vec = 0; vec < all_output_vectors.size(); ++vec) {
      line += '\t' + convert_to_string(all_output_vectors[vec][i], 10);
    }
    // The main data file is always the first entry in the vector of streams
    file_stream << (i + 1) << line << std::endl;
  }
}

void stat_file::time_stamp(std::ofstream &file_stream,
                           std::string const &variables) {
  /**
   * Dates the file and allows the input of a header
   * Input a file stream to write and string of characters to display as
   * headers.
   *
   * @param &stream: Stream that should be timestaped
   * @param variables: A header that can be included underneath the timestamp
   */
  std::chrono::time_point<std::chrono::system_clock> instance;
  instance = std::chrono::system_clock::now();
  std::time_t date_time = std::chrono::system_clock::to_time_t(instance);
  file_stream << "# Created on: " << std::ctime(&date_time);
  file_stream << variables << std::endl;
}

std::string stat_file::file_naming(std::string const &prefix,
                                   size_t const &STEPS, size_t const &N,
                                   double const &DENSITY,
                                   double const &TEMPERATURE,
                                   double const &POWER, double const &A_cst) {
  /**
   * Generates a unique filename for the simulation results to be stored.
   * The method infers from the constructor the number of particles used
   * and the duration of the simulation (steps).
   *
   * @param prefix: File name preceeding the file id
   * @param DENSITY: Fluid density
   * @param TEMPERATURE: Fluid temperature
   * @param POWER: Pair potential power
   * @param A_cst: Softening parameter constant
   *
   * @return: string structured as follows
   *          INPUT_DIR/prefix_step_#_particles_#_rho_#_T_#_n_#_A_#.log
   */

  // Individual streams handling double to string conversions
  std::stringstream A_stream, rho_stream, T_stream;

  rho_stream << std::fixed << std::setprecision(4) << DENSITY;    // 4 decimals
  T_stream << std::fixed << std::setprecision(4) << TEMPERATURE;  // 4 decimals

  _step_to_str = "_step_" + std::to_string(STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();

  // Do not add an A parameter or a potential strength in case they are NAN
  if (isnan(POWER))
    _n_to_str = "";
  else
    _n_to_str = "_n_" + convert_to_string(POWER, 2);

  if (isnan(A_cst))
    _A_to_str = "";
  else {
    A_stream << std::fixed << std::setprecision(5) << A_cst;  // 5 decimals
    _A_to_str = "_A_" + A_stream.str();
  }

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str + _T_to_str +
             _n_to_str + _A_to_str;

  // Explicit defitions
  std::string _FILE_EXT = ".log";

  // Path addition
  return prefix + _FILE_ID + _FILE_EXT;
}

std::string stat_file::convert_to_string(const double &x,
                                         const int &precision) {
  /**
   * Convert doubles to a string with a variable degree of precision.
   *
   * @param &x: Double number to be converted
   * @param &precision: Precision of the double when converted to string
   *
   * @return: string
   */

  std::ostringstream ss;
  ss.str(std::string());  // don't forget to empty the stream
  ss << std::fixed << std::setprecision(precision) << x;

  return ss.str();
}

/**
 * @brief a wrapper for the FileIO::Write2File method
 * 
 * @param output_quantities 
 * @param fstream 
 * @param header 
 */
void stat_file::write_file(std::vector<std::vector<double>> &output_quantities,
                           std::ofstream &fstream, std::string const &header) {
  std::string new_header = "";
  time_stamp(fstream, new_header);
  new_header += header;
  FileIO::Write2File<double>(output_quantities, fstream, "\t", new_header,
                             false);
}
