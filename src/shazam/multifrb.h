#ifndef MULTIFRB_H
#define MULTIFRB_H

#include <cmath>
#include <cstring>
#include <ctime>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/shm.h>
#include <sys/time.h>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "multihdr.h"
#include "utilities.h"

namespace nb = nanobind;
using namespace nb::literals;
using Array = nb::ndarray<nb::numpy, unsigned char, nb::ndim<2>>;

constexpr int FFTSAMPS = 800;
constexpr int FRBMAXBLKS = 12;
constexpr int FRBBUFKEY = 2032;
constexpr long FRBBLKSAMPS = FFTSAMPS * 32;

typedef struct {
  unsigned int active;
  unsigned int status;
  unsigned int empty;
  double pctime;
  double reftime;
  double rectime;
  struct timeval timestamps[FRBMAXBLKS];
  double nanoseconds[FRBMAXBLKS];
  unsigned int flag;
  unsigned int curblk;
  unsigned int currec;
  unsigned int blksize;
  unsigned int nbeams;
  int overflow;
} BeamBufferType;

class MultiFRBSHM {
public:
  MultiFRBSHM()
      : m_header(), m_nf(0), m_nbits(8), m_fh(0.0), m_fl(0.0), m_df(0.0),
        m_bw(0.0), m_dt(0.0), m_mjd(0.0), m_nstokes(1), m_flipped(false),
        m_ra(0.0), m_dec(0.0), m_obsdate(""), m_obstime(""), m_source(""),
        m_beammode(""), m_observer(""), m_gtaccode(""), m_gtactitle(""),
        m_antmaskpol1(0), m_antmaskpol2(0),
        m_antspol1(std::vector<std::string>()),
        m_antspol2(std::vector<std::string>()), m_beamid(0), m_hostid(0),
        m_nbeams(0), m_npcbaselines(0), m_nbeamspernode(0), m_hostname(""),
        m_beamras(std::vector<double>(0.0)),
        m_beamdecs(std::vector<double>(0.0)), m_hdrid(0), m_linked(false),
        m_hdrptr(NULL) {};

  ~MultiFRBSHM() {};

  MultiHeader header() { return m_header; }

  /** Data parameters. **/
  int nf() { return m_nf; };
  double fh() { return m_fh; };
  double fl() { return m_fl; };
  double df() { return m_df; };
  double bw() { return m_bw; };
  double dt() { return m_dt; };
  int nbits() { return m_nbits; };
  int nstokes() { return m_nstokes; };
  bool flipped() { return m_flipped; }

  /** Observation parameters. **/
  double ra() { return m_ra; };
  double dec() { return m_dec; };
  std::string source() { return m_source; };
  std::string obsdate() { return m_obsdate; };
  std::string obstime() { return m_obstime; };
  std::string beammode() { return m_beammode; };
  std::string observer() { return m_observer; };
  std::string gtaccode() { return m_gtaccode; };
  std::string gtactitle() { return m_gtactitle; };
  unsigned int antmaskpol1() { return m_antmaskpol1; };
  unsigned int antmaskpol2() { return m_antmaskpol2; };
  std::vector<std::string> antspol1() { return m_antspol1; };
  std::vector<std::string> antspol2() { return m_antspol2; };

  /** Beam tiling and steering parameters. **/
  int beamid() { return m_beamid; };
  int hostid() { return m_hostid; };
  int nbeams() { return m_nbeams; };
  std::string hostname() { return m_hostname; };
  int npcbaselines() { return m_npcbaselines; };
  int nbeamspernode() { return m_nbeamspernode; };
  std::vector<double> beamras() { return m_beamras; };
  std::vector<double> beamdecs() { return m_beamdecs; };

  /** Shared memory parameters **/
  bool linked() { return m_linked; }
  int maxblks() { return FRBMAXBLKS; }
  int blksamps() { return FRBBLKSAMPS; }

  bool empty() { return m_bufptr->empty; }
  bool status() { return m_bufptr->status; }
  bool active() { return m_bufptr->active; }
  unsigned int curblk() { return m_bufptr->curblk; }
  unsigned int currec() {
    return (m_bufptr->empty) ? m_bufptr->currec
                             : (m_bufptr->currec - 1) % maxblks();
  }

  long blksize() { return blksamps() * m_nf; }
  long size() { return maxblks() * blksize(); }
  double blktime() { return blksamps() * m_dt; }
  double timeofblk(int blk) { return blk * blktime(); }
  double curtime() { return timeofblk(curblk()); }
  double endtime() { return timeofblk(std::ceil(curblk() / maxblks())); }
  double begtime() { return timeofblk(std::floor(curblk() / maxblks())); }

  /** Public methods. **/
  void link();
  void unlink();
  Array getblk(int beam, int blk);
  Array getblks(int beam, int blk0, int blkN);
  Array getslice(int beam, double tbeg, double tend);
  Array getburst(int beam, double t0, double dm, double width);

private:
  /** Shared memory header. **/
  MultiHeader m_header;

  /** Data parameters. **/
  int m_nf;
  int m_nbits;
  double m_fh;
  double m_fl;
  double m_df;
  double m_bw;
  double m_dt;
  double m_mjd;
  int m_nstokes;
  bool m_flipped;

  /** Observation parameters. **/
  double m_ra;
  double m_dec;
  std::string m_obsdate;
  std::string m_obstime;
  std::string m_source;
  std::string m_beammode;
  std::string m_observer;
  std::string m_gtaccode;
  std::string m_gtactitle;
  unsigned int m_antmaskpol1;
  unsigned int m_antmaskpol2;
  std::vector<std::string> m_antspol1;
  std::vector<std::string> m_antspol2;

  /** Beam tiling and steering parameters. **/
  int m_beamid;
  int m_hostid;
  int m_nbeams;
  int m_npcbaselines;
  int m_nbeamspernode;
  std::string m_hostname;
  std::vector<double> m_beamras;
  std::vector<double> m_beamdecs;

  /** Shared memory parameters. **/
  int m_hdrid;
  int m_bufid;
  bool m_linked;
  BeamHeaderType *m_hdrptr;
  BeamBufferType *m_bufptr;
  unsigned char *m_dataptr;

  /** Shared memory pointers. **/
  unsigned char *ptrtobeam(int beam);
  unsigned char *ptrtoblk(int beam, int blk);
  unsigned char *ptrtotime(int beam, double t);
};

void initmultifrb(nb::module_ m);

#endif
