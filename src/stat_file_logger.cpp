#include "stat_file_logger.h"
#include "FileIO.h"

/* Passes the name of the variable to a stream */
#define GET_VAR_NAME(stream, variable) (stream) << #variable

stat_file::stat_file() {}

std::vector<std::ofstream> stat_file::open_files(
    std::vector<std::string> const &file_names) {
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
  /* Write the timestamp and header to the stream */
  time_stamp(file_stream, header);

  /* Delimiter string */
  std::string del = ",";
  /* Find the largest size vector in all our vectors */
  size_t rows = all_output_vectors[0].size();
  for (auto const &i : all_output_vectors)
    if (i.size() > rows) rows = i.size();

  /* This accesses the all_output_vectors with a column major ordering,
     not the best for performance */
  for (size_t i = 0; i < rows; ++i) {
    std::string line = "";
    for (size_t vec = 0; vec < all_output_vectors.size(); ++vec) {
      try {
        line += del + convert_to_string(all_output_vectors.at(vec).at(i), 10);
      }
      /* if the array goes out of bounds then just add 0s */
      catch (const std::out_of_range &e) {
        line += del + "0.0000000000";
      }
    }
    /* The main data file is always the first entry in the vector of streams */
    file_stream << (i + 1) << line << std::endl;
  }
}

void stat_file::time_stamp(std::ofstream &file_stream,
                           std::string const &variables) {
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
  /* Individual streams handling double to string conversions */
  std::stringstream A_stream, rho_stream, T_stream;

  rho_stream << std::fixed << std::setprecision(4) << DENSITY;    // 4 decimals
  T_stream << std::fixed << std::setprecision(4) << TEMPERATURE;  // 4 decimals

  _step_to_str = "_step_" + std::to_string(STEPS);
  _particles_to_str = "_particles_" + std::to_string(N);
  _rho_to_str = "_rho_" + rho_stream.str();
  _T_to_str = "_T_" + T_stream.str();

  /* Do not add an A parameter or a potential strength in case they are NAN */
  if (isnan(POWER))
    _n_to_str = "";
  else
    _n_to_str = "_n_" + convert_to_string(POWER, 2);

  if (isnan(A_cst))
    _A_to_str = "";
  else {
    A_stream << std::fixed << std::setprecision(5) << A_cst; /* 5 decimals */
    _A_to_str = "_A_" + A_stream.str();
  }

  _FILE_ID = _step_to_str + _particles_to_str + _rho_to_str + _T_to_str +
             _n_to_str + _A_to_str;

  /* Explicit defitions */
  std::string _FILE_EXT = ".log";

  /* Path addition */
  return prefix + _FILE_ID + _FILE_EXT;
}

std::string stat_file::convert_to_string(const double &x,
                                         const int &precision) {
  std::ostringstream ss;
  ss.str(std::string()); /* don't forget to empty the stream */
  ss << std::fixed << std::setprecision(precision) << x;

  return ss.str();
}

void stat_file::write_file(std::vector<std::vector<double>> &output_quantities,
                           std::ofstream &fstream, std::string const &header) {
  std::string new_header = "";
  time_stamp(fstream, new_header);
  new_header += header;
  FileIO::Write2File<double>(output_quantities, fstream, ",", new_header,
                             false);
}

void stat_file::write_file_line(std::vector<double> const &output_line,
                                std::ofstream &fstream, int index) {
  /* If the index is no-negative write it in file */
  if (index >= 0) fstream << index << ',';

  std::string line = "";
  for (auto const &i : output_line) line += ',' + convert_to_string(i, 10);
  fstream << line << std::endl;
}
