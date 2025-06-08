#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>
#include <ctime>
#include <memory>
#include <stdexcept>
#include <string>

constexpr double KDM = 1 / 2.41e-4;

double dm2delay(double f, double fref, double dm);
std::string sstrftime(const char *fmt, const std::tm *t);

template <typename... T>
std::string sfmt(const std::string &format, T... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;
  if (size_s <= 0)
    throw std::runtime_error("Could not format string. Exiting...");
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(), buf.get() + size - 1);
}

#endif
