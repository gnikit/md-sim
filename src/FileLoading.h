#pragma once
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#define RESERVE_MEMORY 50000

template <class T>
class FileLoading {
 protected:
 public:
  FileLoading();
  ~FileLoading();
  typedef std::vector<std::vector<T>> vec2d;

  //TODO: make it return a tuple of vectors
  // could be done by either templating a second template variable TUPLE_N
  // ---> TUPLE_N FileLoading<T>::LoadTxt(const std::string & file_name, blah blah blah)
  // or using a vardic template
  //--->http://aherrmann.github.io/programming/2016/02/28/unpacking-tuples-in-cpp14/
  //TODO: Investigate why vec2d not working outside of class def in template LoadTxt
  std::vector<std::vector<T>> LoadTxt(const std::string& file_name,
                                      size_t columns, char comment);
  std::vector<T> LoadSingleCol(const std::string& file_name);
};

template <class T>
FileLoading<T>::FileLoading() {
}

template <class T>
FileLoading<T>::~FileLoading() {
}

template <class T>
std::vector<std::vector<T>> FileLoading<T>::LoadTxt(const std::string& file_name,
                                                    size_t columns, char comment) {
  /*
	Used to read a file with either "space" or "tab" delimated columns.
	The number of columns in the file has to be known and be passed as an argument.

	In the case where a smaller number of columns is supplied, than the number
	contained in the file, then only these columns will be extracted.

	@param file_name: The name/relative path to the file with extension
	@param columns: The total number of columns in the file
	@param comment: Character to be treated as comment. Lines starting with comment
	will be ignored
	@return data[][]: Vector of vectors, with each sub-vector being a read column
	*/

  std::ifstream file(file_name);
  file.exceptions(file.failbit);  // Throws exception
  std::vector<std::vector<T>> data;

  try {
    T temp;
    data.resize(columns);
    std::string line;

    // Reserve space for column vectors
    for (int i = 0; i < data.size(); i++) {
      data.at(i).reserve(RESERVE_MEMORY);
    }

    while (std::getline(file, line)) {
      if (line.size() != 0 && line[0] != comment) {
        std::stringstream ss(line);

        // could include testing for tabs == columns -1
        // otherwise throw range_error or logic_error exception
        // catch that const std::exception &e
        for (int i = 0; i < columns; i++) {
          ss >> temp;
          data[i].push_back(temp);
        }
      }
    }
    file.close();  // just in case it fails to close
  } catch (const std::ios_base::failure& e) {
    if (!file.eof())
      /*if EOF sets eofbit/failbit 0, normally getline sets to 1
			if nothing is extracted from a line. */
      std::cerr << "Caught an ios_base::failure.\n"
                << "Explanatory string: " << e.what() << '\n'
                << "Error code: " << e.code() << '\n';
  }
  return data;
}

template <class T>
std::vector<T> FileLoading<T>::LoadSingleCol(const std::string& file_name) {
  /*
	Reads a file that is structured with data in a column
	*/
  std::vector<T> data;
  data.reserve(RESERVE_MEMORY);  // increases performance for large files
  std::ifstream file(file_name);
  file.exceptions(file.failbit);  // Throws exception

  try {
    std::copy(std::istream_iterator<T>(file),
              std::istream_iterator<T>(), std::back_inserter(data));
  }

  catch (const std::ios_base::failure& e) {
    if (!file.eof())
      /*if EOF sets eofbit/failbit 0, normally getline sets to 1
			if nothing is extracted from a line. */
      std::cerr << "Caught an ios_base::failure.\n"
                << "Explanatory string: " << e.what() << '\n'
                << "Error code: " << e.code() << '\n';
  }
  return data;
}
