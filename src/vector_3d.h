/**
 * @file vector_3d<T>.h
 * @author Ioannis Nikiteas
 * @brief This is just a convinient little template to store 3d data
 * @version 0.1
 * @date 2020-05-07
 *
 * @copyright Copyright (c) 2020
 *
 */

#pragma once
#include <iostream>
#include <vector>

using namespace std;

/**
 * @brief Class which holds data for x, y, z data.
 * The default (implicit) operator= and copy constructors are being used.
 *
 */
template <typename T>
class vector_3d {
 private:
  bool xyz_same_size = true;

 public:
  vector<T> x, y, z;

  /******************************* CONSTRUCTORS *******************************/

  vector_3d<T>();
  vector_3d<T>(vector<T> X, vector<T> Y, vector<T> Z);

  /*************************** OPERATOR OVERLOADING ***************************/

  vector_3d<T>& operator+=(vector_3d<T> const& rhs);
  vector_3d<T>& operator+=(T const& rhs);

  vector_3d<T>& operator-=(vector_3d<T> const& rhs);
  vector_3d<T>& operator-=(T const& rhs);

  vector_3d<T>& operator*=(vector_3d<T> const& rhs);
  vector_3d<T>& operator*=(T const& rhs);

  vector_3d<T>& operator/=(vector_3d<T> const& rhs);
  vector_3d<T>& operator/=(T const& rhs);

  /*************************** METHODS OVERLOADING ****************************/

  /**
   * @brief Returns x.size() assuming x.size() == y.size() == z.size()
   *
   * @return size_t vector size
   */
  size_t size();

  void resize(size_t const& sz);
  void resize(size_t const& sz, T val);
  void resize(size_t const& szx, size_t const& szy, size_t const& szz);

  void reserve(size_t const& sz);
  void reserve(size_t const& szx, size_t const& szy, size_t const& szz);

  /********************************* METHODS **********************************/

  /**
   * @brief Calculates the magnitude at each point of the vector_3d<T>
   * magnitude = sqrt(x^2 + y^2 + z^2)
   *
   * @return std::vector<T> Magnitude
   */
  std::vector<T> magnitude();
};

template <typename T>
vector_3d<T> operator+(vector_3d<T> const& lhs, vector_3d<T> const& rhs);
template <typename T>
vector_3d<T> operator+(vector_3d<T> const& lhs, T const& rhs);

template <typename T>
vector_3d<T> operator-(vector_3d<T> const& lhs, vector_3d<T> const& rhs);
template <typename T>
vector_3d<T> operator-(vector_3d<T> const& lhs, T const& rhs);

template <typename T>
vector_3d<T> operator*(vector_3d<T> const& lhs, vector_3d<T> const& rhs);
template <typename T>
vector_3d<T> operator*(vector_3d<T> const& lhs, T const& rhs);

template <typename T>
vector_3d<T> operator/(vector_3d<T> const& lhs, vector_3d<T> const& rhs);
template <typename T>
vector_3d<T> operator/(vector_3d<T> const& lhs, T const& rhs);

template <typename T>
std::ostream& operator<<(std::ostream& out, const vector_3d<T>& v);

/*******************************************************************************
 * "vector_3d.cpp" code follows
 * ****************************************************************************/

/******************************** CONSTRUCTORS ********************************/
template <typename T>
vector_3d<T>::vector_3d() {}

template <typename T>
vector_3d<T>::vector_3d(vector<T> X, vector<T> Y, vector<T> Z)
    : x{X}, y{Y}, z{Z} {
  if (!(x.size() == y.size() && x.size() == z.size())) xyz_same_size = false;
}

/**************************** OPERATOR OVERLOADING ****************************/

template <typename T>
vector_3d<T>& vector_3d<T>::operator+=(vector_3d<T> const& rhs) {
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
template <typename T>
vector_3d<T>& vector_3d<T>::operator+=(T const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] += rhs;
      this->y[i] += rhs;
      this->z[i] += rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, T" << std::endl;
    exit(-1);
  }

  return *this;
}

template <typename T>
vector_3d<T>& vector_3d<T>::operator-=(vector_3d<T> const& rhs) {
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
template <typename T>
vector_3d<T>& vector_3d<T>::operator-=(T const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] -= rhs;
      this->y[i] -= rhs;
      this->z[i] -= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, T" << std::endl;
    exit(-1);
  }

  return *this;
}

template <typename T>
vector_3d<T>& vector_3d<T>::operator*=(vector_3d<T> const& rhs) {
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
template <typename T>
vector_3d<T>& vector_3d<T>::operator*=(T const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] *= rhs;
      this->y[i] *= rhs;
      this->z[i] *= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, T" << std::endl;
    exit(-1);
  }

  return *this;
}

template <typename T>
vector_3d<T>& vector_3d<T>::operator/=(vector_3d<T> const& rhs) {
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
template <typename T>
vector_3d<T>& vector_3d<T>::operator/=(T const& rhs) {
  if (this->xyz_same_size) {
    for (size_t i = 0; i < this->size(); ++i) {
      this->x[i] /= rhs;
      this->y[i] /= rhs;
      this->z[i] /= rhs;
    }
  } else {
    std::cerr << "Need to overload '=' for vector, T" << std::endl;
    exit(-1);
  }

  return *this;
}

/**************************** METHODS OVERLOADING *****************************/

template <typename T>
size_t vector_3d<T>::size() {
  if (!(x.size() == y.size() && x.size() == z.size())) abort();

  return x.size();
}

template <typename T>
void vector_3d<T>::resize(size_t const& sz) {
  x.resize(sz);
  y.resize(sz);
  z.resize(sz);
}

template <typename T>
void vector_3d<T>::resize(size_t const& sz, T val) {
  x.resize(sz, val);
  y.resize(sz, val);
  z.resize(sz, val);
}

template <typename T>
void vector_3d<T>::resize(size_t const& szx, size_t const& szy,
                          size_t const& szz) {
  x.resize(szx);
  y.resize(szy);
  z.resize(szy);
}

template <typename T>
void vector_3d<T>::reserve(size_t const& sz) {
  x.reserve(sz);
  y.reserve(sz);
  z.reserve(sz);
}

template <typename T>
void vector_3d<T>::reserve(size_t const& szx, size_t const& szy,
                           size_t const& szz) {
  x.reserve(szx);
  y.reserve(szy);
  z.reserve(szy);
}

/********************************** METHODS ***********************************/

template <typename T>
std::vector<T> vector_3d<T>::magnitude() {
  try {
    std::vector<T> magn(x.size());

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

/************************ GLOBAL OPERATOR OVERLOADS ***************************/

template <typename T>
vector_3d<T> operator+(vector_3d<T> const& lhs, vector_3d<T> const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d<T> tmp(lhs);
  tmp += rhs;
  return tmp;
}
template <typename T>
vector_3d<T> operator+(vector_3d<T> const& lhs, T const& rhs) {
  vector_3d<T> tmp(lhs);
  tmp += rhs;
  return tmp;
}

template <typename T>
vector_3d<T> operator-(vector_3d<T> const& lhs, vector_3d<T> const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d<T> tmp(lhs);
  tmp -= rhs;
  return tmp;
}
template <typename T>
vector_3d<T> operator-(vector_3d<T> const& lhs, T const& rhs) {
  vector_3d<T> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

template <typename T>
vector_3d<T> operator*(vector_3d<T> const& lhs, vector_3d<T> const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d<T> tmp(lhs);
  tmp *= rhs;
  return tmp;
}
template <typename T>
vector_3d<T> operator*(vector_3d<T> const& lhs, T const& rhs) {
  vector_3d<T> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

template <typename T>
vector_3d<T> operator/(vector_3d<T> const& lhs, vector_3d<T> const& rhs) {
  /* Call the IMPLICITLY generated copy constructor of the compiler */
  vector_3d<T> tmp(lhs);
  tmp /= rhs;
  return tmp;
}
template <typename T>
vector_3d<T> operator/(vector_3d<T> const& lhs, T const& rhs) {
  vector_3d<T> tmp(lhs);
  tmp /= rhs;
  return tmp;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, vector_3d<T> const& v) {
  for (size_t i = 0; i < v.x.size(); ++i)
    out << v.x[i] << ", " << v.y[i] << ", " << v.z[i] << std::endl;

  return out;
}
