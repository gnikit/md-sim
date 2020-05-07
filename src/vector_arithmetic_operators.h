/**
 * @file vector_arithmetic_operators.h
 * @author Ioannis Nikiteas
 * @brief Simple arithmetic overloads for std::vector. If this
 * https://github.com/cplusplus/papers/issues/169 ever comes to fruition
 * and joins the standard this template will be unnecessary and will have to be
 * deleted
 * @date 2020-05-07
 *
 */

#pragma once
#include <iostream>  /* std::cerr */
#include <stdexcept> /* std::out_of_range */
#include <vector>    /* std::vector */

template <typename T>
std::vector<T> operator+(std::vector<T> const& lhs, std::vector<T> const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] += rhs[i];
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }

  return tmp;
}

template <typename T>
T operator+(std::vector<T> const& lhs, T const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] += rhs;
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}

template <typename T>
std::vector<T> operator-(std::vector<T> const& lhs, std::vector<T> const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] -= rhs[i];
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }

  return tmp;
}

template <typename T>
T operator-(std::vector<T> const& lhs, T const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] -= rhs;
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}

template <typename T>
std::vector<T> operator*(std::vector<T> const& lhs, std::vector<T> const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] *= rhs[i];
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}

template <typename T>
T operator*(std::vector<T> const& lhs, T const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] *= rhs;
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}

template <typename T>
std::vector<T> operator/(std::vector<T> const& lhs, std::vector<T> const& rhs) {
  std::vector<T> tmp(lhs);

  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] /= rhs[i];
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}
template <typename T>
T operator/(std::vector<T> const& lhs, T const& rhs) {
  std::vector<T> tmp(lhs);
  try {
    for (size_t i = 0; i < tmp.size(); ++i) tmp[i] /= rhs;
  } catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    exit(-1);
  }
  return tmp;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  out << "[";
  size_t i = 0;
  for (i = 0; i < v.size() - 1; ++i) out << v[i] << ", ";

  out << v[i] << "]" << std::endl;

  return out;
}