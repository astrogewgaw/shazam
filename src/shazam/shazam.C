#include <nanobind/nanobind.h>

#include "multifrb.h"
#include "shazam.h"

NB_MODULE(_internals, m) { init_multifrb(m); }
