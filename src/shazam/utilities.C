#include "utilities.h"

double dm2delay(double f, double fref, double dm) {
  return KDM * dm * (std::pow(f, -2) - std::pow(fref, -2));
}
