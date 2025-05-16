#include "multifrb.h"
#include "nanobind/nanobind.h"

class MultiTELSHM {};
class MultiFRBSHM {};

void init_multifrb(nb::module_ m) {
  nb::class_<MultiFRBSHM>(m, "MultiFRBSHM").def(nb::init<>());
}
