#include "data_structures.h"

/******************************** CONSTRUCTORS ********************************/

vector_3d::vector_3d() {}

vector_3d::vector_3d(vector<double> X, vector<double> Y, vector<double> Z)
    : x{X}, y{Y}, z{Z} {
  if (!(x.size() == y.size() && x.size() == z.size())) xyz_same_size = false;
}

/**************************** OPERATOR OVERLOADING ****************************/

vector_3d& vector_3d::operator+=(vector_3d const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] += rhs.x[i];
      this->y[i] += rhs.y[i];
      this->z[i] += rhs.z[i];
    }
  } else {
    this->x = this->x + rhs.x;
    this->y = this->y + rhs.y;
    this->z = this->z + rhs.z;
  }

  return *this;
}
vector_3d& vector_3d::operator+=(double const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] += rhs;
      this->y[i] += rhs;
      this->z[i] += rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, double" << std::endl;
    exit(-1);
  }

  return *this;
}

vector_3d& vector_3d::operator-=(vector_3d const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] -= rhs.x[i];
      this->y[i] -= rhs.y[i];
      this->z[i] -= rhs.z[i];
    }
  } else {
    this->x = this->x - rhs.x;
    this->y = this->y - rhs.y;
    this->z = this->z - rhs.z;
  }

  return *this;
}
vector_3d& vector_3d::operator-=(double const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] -= rhs;
      this->y[i] -= rhs;
      this->z[i] -= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, double" << std::endl;
    exit(-1);
  }

  return *this;
}

vector_3d& vector_3d::operator*=(vector_3d const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] *= rhs.x[i];
      this->y[i] *= rhs.y[i];
      this->z[i] *= rhs.z[i];
    }
  } else {
    this->x = this->x * rhs.x;
    this->y = this->y * rhs.y;
    this->z = this->z * rhs.z;
  }

  return *this;
}
vector_3d& vector_3d::operator*=(double const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] *= rhs;
      this->y[i] *= rhs;
      this->z[i] *= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, double" << std::endl;
    exit(-1);
  }

  return *this;
}

vector_3d& vector_3d::operator/=(vector_3d const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] /= rhs.x[i];
      this->y[i] /= rhs.y[i];
      this->z[i] /= rhs.z[i];
    }
  } else {
    this->x = this->x / rhs.x;
    this->y = this->y / rhs.y;
    this->z = this->z / rhs.z;
  }

  return *this;
}
vector_3d& vector_3d::operator/=(double const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] /= rhs;
      this->y[i] /= rhs;
      this->z[i] /= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, double" << std::endl;
    exit(-1);
  }

  return *this;
}

vector_3d& operator+(vector_3d const& lhs, vector_3d const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d tmp(lhs);
  tmp += rhs;
  return tmp;
}
vector_3d& operator+(vector_3d const& lhs, double const& rhs) {
  vector_3d tmp(lhs);
  tmp += rhs;
  return tmp;
}

vector_3d& operator-(vector_3d const& lhs, vector_3d const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d tmp(lhs);
  tmp -= rhs;
  return tmp;
}
vector_3d& operator-(vector_3d const& lhs, double const& rhs) {
  vector_3d tmp(lhs);
  tmp -= rhs;
  return tmp;
}

vector_3d& operator*(vector_3d const& lhs, vector_3d const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d tmp(lhs);
  tmp *= rhs;
  return tmp;
}
vector_3d& operator*(vector_3d const& lhs, double const& rhs) {
  vector_3d tmp(lhs);
  tmp *= rhs;
  return tmp;
}

vector_3d& operator/(vector_3d const& lhs, vector_3d const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d tmp(lhs);
  tmp /= rhs;
  return tmp;
}
vector_3d& operator/(vector_3d const& lhs, double const& rhs) {
  vector_3d tmp(lhs);
  tmp /= rhs;
  return tmp;
}

std::ostream& operator<<(std::ostream& out, vector_3d const& v) {
  for (size_t i = 0; i < v.x.size(); ++i)
    out << v.x[i] << ", " << v.y[i] << ", " << v.z[i] << std::endl;

  return out;
}

/**************************** METHODS OVERLOADING *****************************/

size_t vector_3d::size() {
  if (!(x.size() == y.size() && x.size() == z.size())) abort();

  return x.size();
}

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

/********************************** METHODS ***********************************/

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