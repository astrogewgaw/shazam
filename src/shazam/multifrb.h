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

namespace nb = nanobind;
using namespace nb::literals;
using Array = nb::ndarray<nb::numpy, unsigned char, nb::ndim<2>>;

constexpr int FRBHDRKEY = 2031;
constexpr int FRBBUFKEY = 2032;

constexpr int FFTSAMPS = 800;
constexpr int FRBMAXBLKS = 12;
constexpr long FRBBLKSAMPS = FFTSAMPS * 32;

constexpr char ANTENNAS[30][4] = {
    "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C08", "C09", "C10",
    "C11", "C12", "C13", "C14", "E02", "E03", "E04", "E05", "E06", "S01",
    "S02", "S03", "S04", "S06", "W01", "W02", "W03", "W04", "W05", "W06"};

constexpr char BEAMTYPES[7][6] = {"IA",  "PA",   "VLT", "PC",
                                  "CDP", "PASV", "MISC"};

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

template <typename... Args>
std::string safe_fmt(const std::string &format, Args... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;
  if (size_s <= 0)
    throw std::runtime_error("Could not format string. Exiting...");
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(), buf.get() + size - 1);
}

inline std::string safe_strftime(const char *fmt, const std::tm *t) {
  std::size_t len = sizeof(fmt);
  auto buff = std::make_unique<char[]>(len);
  while (std::strftime(buff.get(), len, fmt, t) == 0) {
    len *= 2;
    buff = std::make_unique<char[]>(len);
  }
  return std::string{buff.get()};
}

class MultiFRBSHM {
public:
  MultiFRBSHM() { link(); };
  ~MultiFRBSHM() { unlink(); };

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
  int maxblks() { return FRBMAXBLKS; }
  int blksamps() { return FRBBLKSAMPS; }
  unsigned int currec() {
    m_currec = (m_bufptr->empty) ? m_bufptr->currec
                                 : (m_bufptr->currec - 1) % maxblks();
    return m_currec;
  }
  unsigned int curblk() {
    m_curblk = m_bufptr->curblk;
    return m_curblk;
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
  Array getslice(int beam, double tbeg, double tend);
  Array getburst(int beam, double t0, double dm, double width);

private:
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
  unsigned int m_currec;
  unsigned int m_curblk;
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
