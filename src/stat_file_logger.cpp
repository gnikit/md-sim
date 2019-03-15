#include "stat_file_logger.h"

// Passes the name of the variable to a stream
#define GET_VAR_NAME(stream, variable) (stream) << #variable

stat_file::stat_file() {}

void stat_file::open_files(std::string &data, std::string &rdf,
                           std::string &pos) {
  /*
   * Open/Create if file does not exist.
   * Overwrites existing data.
   */

  DATA.open(data, std::ios::out | std::ios::trunc);
  RDF.open(rdf, std::ios::out | std::ios::trunc);
  POS.open(pos, std::ios::out | std::ios::trunc);
}

void stat_file::write_data_file(
    size_t STEPS, std::vector<double> &density,
    std::vector<double> &temperature, std::vector<double> &u_en,
    std::vector<double> &k_en, std::vector<double> &pc, std::vector<double> &pk,
    std::vector<double> &msd, std::vector<double> &Cr, std::vector<double> &sfx,
    std::vector<double> &sfy, std::vector<double> &sfz) {
  /*
   * Writes values of parameters to Data file.
   */
  for (size_t i = 0; i < STEPS; ++i) {
    DATA << (i + 1) << '\t' << density[i] << '\t' << temperature[i] << '\t'
         << u_en[i] << '\t' << k_en[i] << '\t' << pc[i] << '\t' << pk[i] << '\t'
         << msd[i] << '\t' << Cr[i] << sfx[i] << '\t' << sfy[i] << '\t'
         << sfz[i] << std::endl;
  }
}

void stat_file::time_stamp(std::ofstream &stream, std::string variables) {
  /*
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
  stream << "# Created on: " << std::ctime(&date_time);
  stream << variables << std::endl;
}

std::string stat_file::file_naming(std::string prefix, size_t &STEPS, size_t &N,
                                   double &DENSITY, double &TEMPERATURE,
                                   double &POWER, double &A_cst) {
  /*
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
   *          INPUT_DIR/prefix_step_#_particles_#_rho_#_T_#_n_#_A_#.txt
   */

  // Individual streams handling double to string conversions
  std::stringstream A_stream, rho_stream, T_stream;

  rho_stream << std::fixed << std::setprecision(4) << DENSITY;    // 4 decimals
  T_stream << std::fixed << std::setprecision(4) << TEMPERATURE;  // 4 decimals
  A_stream << std::fixed << std::setprecision(5) << A_cst;        // 5 decimals

  _step_to_str = "_step_" + std::to_string(STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();
  _n_to_str = "_n_" + convert_to_string(POWER, 2);
  _A_to_str = "_A_" + A_stream.str();

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str + _T_to_str +
             _n_to_str + _A_to_str;

  // Explicit defitions
  std::string _FILE_EXT = ".txt";

  // Path addition
  return prefix + _FILE_ID + _FILE_EXT;
}
std::string stat_file::convert_to_string(const double &x,
                                         const int &precision) {
  /*
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
