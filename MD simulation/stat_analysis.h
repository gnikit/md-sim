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

private:

  vec1d TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6;
  vec1d _A_list;
  std::string _path;

  long double _temp0 = 0;
  long double _temp1 = 0;
  long double _temp2 = 0;
  long double _temp3 = 0;
  long double _temp4 = 0;
  long double _temp5 = 0;
  long double _temp6 = 0;

public:
  Stat_Analysis(std::string PATH, vec1d A_LIST);
  ~Stat_Analysis();

  void ReadFromFile(vec1d &x, const std::string &file_name);
  void ReadFromFile(vec1d &x, vec1d &y, vec1d &z,
                    vec1d &Xx, vec1d &Yy, vec1d &Zz, vec1d &w,
                    const std::string &file_name);
  void Mean(vec1d &x, vec1d &y, vec1d &z,
            vec1d &Xx, vec1d &Yy, vec1d &Zz, vec1d &w);
  void StaticDataProcessing(size_t n);
  // TODO: make a RunMe method and add arguments to the existing methods for increased reusability

};