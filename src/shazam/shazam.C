#include "multifrb.h"
#include "multihdr.h"

NB_MODULE(internals, m) {
  initmultihdr(m);
  initmultifrb(m);
}
