#include "multifrb.h"
#include "multihdr.h"
#include "multitel.h"

NB_MODULE(internals, m) {
  initmultihdr(m);
  initmultitel(m);
  initmultifrb(m);
}
