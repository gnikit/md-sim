#include "data_structures.h"

void vector_3d::resize(size_t const& sz) {
  x.resize(sz);
  y.resize(sz);
  z.resize(sz);
}

void vector_3d::resize(size_t const& sz, double val) {
  x.resize(sz, val);
  y.resize(sz, val);
  z.resize(sz, val);
}

void vector_3d::resize(size_t const& szx, size_t const& szy,
                       size_t const& szz) {
  x.resize(szx);
  y.resize(szy);
  z.resize(szy);
}

void vector_3d::reserve(size_t const& sz) {
  x.reserve(sz);
  y.reserve(sz);
  z.reserve(sz);
}

void vector_3d::reserve(size_t const& szx, size_t const& szy,
                        size_t const& szz) {
  x.reserve(szx);
  y.reserve(szy);
  z.reserve(szy);
}

std::vector<double> vector_3d::magnitude() {
  try {
    std::vector<double> magn(x.size());

    for (size_t i = 0; i < x.size(); ++i)
      // todo: maybe use std::hypot
      magn[i] = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

    return magn;

  } catch (const std::out_of_range& e) {
    std::cerr << "Out of range exception in magnitude calculation.\n"
              << "This probably means that x,y,z are not te same size!\n"
              << "Error message: " << e.what() << std::endl;
    exit(-1);
  }
}