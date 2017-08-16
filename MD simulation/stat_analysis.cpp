#include "stat_analysis.h"

Stat_Analysis::Stat_Analysis(std::string PATH, vec1d A_LIST) {
  // Do everything in the constructor
  // No need for a Run function
  _path = PATH;
  _A_list = A_LIST;
}

Stat_Analysis::~Stat_Analysis() {
}


void Stat_Analysis::ReadFromFile(vec1d &x, const std::string &file_name) {
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

void Stat_Analysis::ReadFromFile(vec1d &x, vec1d &y, vec1d &z,
  vec1d &Xx, vec1d &Yy, vec1d &Zz, vec1d &w,
  const std::string &file_name) {
  /*
  Reads from a stream that already exists for a file that is already placed in the
  directory and appends the data into a 1D vector.
  */
  std::ifstream read_file(file_name);
  try {
    if (read_file.is_open() == false) {
      throw "File does not exist!";
    }

    //assert(read_file.is_open()); // Evaluates false if file does not exist
    long double a, b, c, d, e, f, g;
    std::string line;
    while (!read_file.eof()) {
      std::getline(read_file, line);

      if (line.length() == 0 || line[0] == '#') {
        std::cout << "IGNORE LINE" << std::endl;
      }
      else {
        std::stringstream ss(line);
        ss >> a;
        ss >> b;
        ss >> c;
        ss >> d;
        ss >> e;
        //ss >> f;
        //ss >> g;
        //TODO: convert string STEP_NUM to int and resize vec1d to the appropriate size
        // add an index ++i that occupies vec1ds
        x.push_back(a);
        y.push_back(b);
        z.push_back(c);
        Xx.push_back(d);
        Yy.push_back(e);
        //Zz.push_back(f);
        //w.push_back(g);
      }
    }
  }
  catch (std::string e) {
    std::cout << e << std::endl;
  }
  read_file.close();

}

void Stat_Analysis::Mean(vec1d &x, vec1d &y, vec1d &z,
  vec1d &Xx, vec1d &Yy, vec1d &Zz, vec1d &w) {
  size_t n = 0;

  for (size_t i = 0; i < x.size(); i++) {
    _temp0 += x[i];
    _temp1 += y[i];
    _temp2 += z[i];
    _temp3 += Xx[i];
    _temp4 += Yy[i];
    //_temp5 += Zz[i];
    //_temp6 += w[i];
    n++;
  }
  _temp0 /= n; // this is the mean
  _temp1 /= n;
  _temp2 /= n;
  _temp3 /= n;
  _temp4 /= n;
  //_temp5 /= n;
  //_temp6 /= n;

}

void Stat_Analysis::StaticDataProcessing(size_t n)
/*
Used to average the energies and pressures for all
values of A for a specific number n of the potential
strength
*/
{
  vec1d temp;

  std::ofstream data;
  std::string file_name, path, sep, power, a, txt;
  path = _path;
  sep = "~";
  txt = ".txt";
  // path = "\0";
  power = std::to_string(n);
  std::string name = path + "AVGdata" + power + txt;
  data.open(name, std::ios::out | std::ios::trunc);

  // data << "Power:\t" << n << std::endl;
  data << "# A:\tT\tK:\tU:\tETot:\tPc:" << std::endl;

  for (size_t i = 0; i < _A_list.size(); i++) {
    file_name.clear();
    TEMP0.clear();
    TEMP1.clear();
    TEMP2.clear();
    TEMP3.clear();
    TEMP4.clear();
    TEMP5.clear();
    TEMP6.clear();
    file_name = "Data";

    // Converts A to 2 decimal long double
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << _A_list.at(i);

    a = stream.str();
    file_name = path + file_name + power + sep + a + txt;



    ReadFromFile(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, file_name);


    Mean(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6);
    data.precision(5);
    data << _A_list.at(i) << '\t' << _temp0 << '\t' << _temp1 << '\t' <<
      _temp2 << '\t' << _temp3 << '\t' << _temp4 << std::endl;
    // consider excluding usless vectors like PK, PTOT
    // Debug to see if the smethod for file reading still works

  }
}