#include "utilities.h"

double dm2delay(double f, double fref, double dm) {
  return KDM * dm * (std::pow(f, -2) - std::pow(fref, -2));
}

std::string sstrftime(const char *fmt, const std::tm *t) {
  std::size_t len = sizeof(fmt);
  auto buff = std::make_unique<char[]>(len);
  while (std::strftime(buff.get(), len, fmt, t) == 0) {
    len *= 2;
    buff = std::make_unique<char[]>(len);
  }
  return std::string{buff.get()};
}
