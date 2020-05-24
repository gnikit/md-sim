#pragma once
#include <math.h> /* isnan */
#include <chrono>
#include <fstream>
#include <iomanip> /* setprecision */
#include <sstream> /* stringstream */
#include <string>
#include <vector>

class stat_file {
 protected:
  const std::string __file_ext = ".csv";
  std::string _dir, __file_id;
  std::string __step2str, __particles2str, __rho2str, __T2str, __n2str, __a2str;

 public:
  std::vector<std::string> file_names;

  /* Initialise them before doing any IO operations */
  stat_file();

  /**
   * @brief creates an array of output streams given a vector of file names
   *
   * @param file_names: vector containing all the filenames
   * @return std::vector<std::ofstream>: vector containing all the output
   * streams
   */
  static std::vector<std::ofstream> open_files(
      std::vector<std::string> const &file_names);

  /**
   * @brief Writes all the output vectors to a file using an existing
   * and open stream. An optional header for the file can be included.
   *
   * @param file_stream: A open out filestream
   * @param header: header string for the file
   * @param all_output_vectors: a vector of vectors containing all the data
   */
  static void write_data_file(
      std::ofstream &file_stream, std::string const &header,
      std::vector<std::vector<double>> const &all_output_vectors);

  /**
   * @brief
   * Dates the file and allows the input of a header
   * Input a file stream to write and string of characters to display as
   * headers.
   *
   * @param &stream: Stream that should be timestaped
   * @param variables: A header that can be included underneath the timestamp
   */
  static void time_stamp(std::ofstream &file_stream,
                         std::string const &variables);

  /**
   * @brief
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
   *          INPUT_DIR/prefix_step_#_particles_#_rho_#_T_#_n_#_A_#.csv
   */
  std::string file_naming(std::string const &prefix, size_t const &STEPS,
                          size_t const &N, double const &DENSITY,
                          double const &TEMPERATURE, double const &POWER,
                          double const &A_cst);

  /**
   *  @brief Convert doubles to a string with a variable degree of precision.
   *
   * @param &x: Double number to be converted
   * @param &precision: Precision of the double when converted to string
   *
   * @return: string
   */
  static std::string convert_to_string(const double &x, const int &precision);

  /**
   * @brief a wrapper for the FileIO::Write2File method
   *
   * @param output_quantities: a 2D vector
   * @param fstream: the file stream
   * @param header: an optional header for the file
   */
  void write_file(std::vector<std::vector<double>> &output_quantities,
                  std::ofstream &fstream, std::string const &header,
                  size_t format = 1);

  void write_file(std::vector<std::vector<double> *> &output_quantities,
                  std::ofstream &fstream, std::string const &header,
                  size_t format = 1);

  /**
   * @brief Write the output_line vector as a line to fstream (ends the line
   * with and endl character). Also checks if fstream has been opened and if not
   * it opens it.
   *
   * @param output_line: Double vector containing the line data
   * @param fstream: Output file stream
   * @param index: (optional) The row index for the line. If negative ignore
   */
  void write_file_line(std::vector<double> const &output_line,
                       std::ofstream &fstream, int index = -1);
};