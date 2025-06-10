#ifndef MULTITEL_H
#define MULTITEL_H

#include <sys/time.h>

#include "multihdr.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAXBLKS 8

enum {
  MaxPols = 4,
  MaxSamples = (16 * MAXBLKS),
  Channels = CHANNEL,
  WordSize = sizeof(short),
  WordMask = 0xffff,
  TotalWords = MaxPols * MaxSamples * Channels,
  RecSize = TotalWords * WordSize,
};

enum {
  UnInitialized = 1 << 0,
  Marked = 1 << 1,
  GoodData = 1 << 2,
  AcqOver = 1 << 3,
  GPSpresent = 1 << 8,
  BlockErr = 1 << 13,
  TimeErr = 1 << 14,
  SimData = 1 << 15
};

enum {
  MaxGPS = 16,
  MaxBLK = 64,
  MaxRecs = 8,
  ExtraWords = 32,
  ExtraBuf = ExtraWords * WordSize
};

enum {
  PC2IST,
  PCref,
  ISTref,
  TimeParams,
  PC2BLK = PC2IST,
  BLKref = ISTref,
  SeqParams = TimeParams
};

enum { SetCount = 100, SetTime, GetTime };

typedef struct {
  unsigned short flag;
  unsigned short count;
  int seq;
  long tv[2];
} TimeType;

enum { TimeSz = sizeof(TimeType) };

typedef struct {
  unsigned data_flag;
  unsigned dbuf_seq;
  unsigned time_err;
  unsigned blocks;
  double prev_time;
  TimeType block_time;
} RackInfoType;

typedef struct {
  unsigned rec_flag;
  unsigned rec_seq;
  unsigned beg_off;
  unsigned short *begp;
  double pc_time;
  double rec_time;
  struct timeval timestamp_gps;
  double blk_nano;
  int AcqSeqNo;
  struct timeval timestamp_pc;
} RecType;

typedef struct {
  unsigned acq_flag;
  unsigned mark_num;
  unsigned dummy;
  void *shmp;
  long ref_time;
  long tzoff;
  double blk_time;
  double pc2ist[TimeParams];
  double pc2blk[SeqParams];
  double gps_val[MaxGPS];
  double blk_val[MaxBLK];
  unsigned gps_seq[MaxGPS];
  unsigned blk_seq[MaxBLK];
  int gps_ind, blk_ind;
  RecType rec[MaxRecs + 1];
  int rec_ind;
  double RecStartTime[MaxRecs + 1];
  int marker_offsets[3][5];
} GlobalInfoType;

enum {
  ShmKey = 1034,
  PageSize = 4096,
  ShmDataOff = ((sizeof(GlobalInfoType) + ExtraBuf) / PageSize + 1) * PageSize,
  ShmDataSize =
      (((MaxRecs + 1) * RecSize + ExtraBuf) / PageSize + 1) * PageSize,
  ShmSize = ShmDataOff + ShmDataSize
};

#define FFT_CYCLE CHANNEL * 2
#define BLK_INTEG 8192
#define BASE_CLK 400.0e6

#ifdef __cplusplus
}
#endif

#include <cstring>
#include <ctime>
#include <stdexcept>
#include <string>
#include <sys/shm.h>
#include <sys/time.h>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using Array = nb::ndarray<nb::numpy, unsigned char, nb::ndim<2>>;

class MultiTELSHM {
public:
  MultiTELSHM()
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

  ~MultiTELSHM() {};

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
  int maxblks() { return MaxRecs; }
  int blksamps() { return 32 * 25; }

  unsigned int curblk() { return m_bufptr->rec[currec()].rec_seq; }
  unsigned int currec() { return (m_bufptr->rec_ind - 1) % maxblks(); }

  long blksize() { return blksamps() * nf(); }
  long size() { return maxblks() * blksize(); }
  double blktime() { return blksamps() * dt(); }
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
  GlobalInfoType *m_bufptr;
  unsigned char *m_dataptr;

  /** Shared memory pointers. **/
  unsigned char *ptrtobeam(int beam);
  unsigned char *ptrtoblk(int beam, int blk);
  unsigned char *ptrtotime(int beam, double t);
};

void initmultitel(nb::module_ m);

#endif
