#ifndef SHAZAM_FRB_H
#define SHAZAM_FRB_H

#include <sys/shm.h>

#include <chrono>
#include <cmath>
#include <string>
#include <tuple>

#include "hdr.h"

namespace shazam {
  constexpr int FRBFFTSAMPS = 800;
  constexpr int FRBMAXBLKS = 12;
  constexpr int FRBBUFKEY = 2032;
  constexpr long FRBBLKSAMPS = FRBFFTSAMPS * 32;

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

  class FRBRing {
  public:
    FRBRing()
        : m_hdr(),
          m_hdrid(0),
          m_bufid(0),
          m_mode(READ),
          m_hdrptr(NULL),
          m_bufptr(NULL),
          m_opened(false),
          m_dataptr(NULL) {}

    ~FRBRing() {}

    MODE mode() { return m_mode; }
    Header hdr() { return m_hdr; }
    bool opened() { return m_opened; }

    /** Data parameters. **/
    int nf() { return m_hdr.m_nf; }
    double fh() { return m_hdr.m_fh; }
    double fl() { return m_hdr.m_fl; }
    double df() { return m_hdr.m_df; }
    double bw() { return m_hdr.m_bw; }
    double dt() { return m_hdr.m_dt; }
    int nbits() { return m_hdr.m_nbits; }
    int nstokes() { return m_hdr.m_nstokes; }
    bool flipped() { return m_hdr.m_flipped; }

    /** Observation parameters. **/
    double ra() { return m_hdr.m_ra; }
    double dec() { return m_hdr.m_dec; }
    std::string source() { return m_hdr.m_source; }
    std::string beammode() { return m_hdr.m_beammode; }
    std::string observer() { return m_hdr.m_observer; }
    std::string gtaccode() { return m_hdr.m_gtaccode; }
    std::string gtactitle() { return m_hdr.m_gtactitle; }
    unsigned int antmaskpol1() { return m_hdr.m_antmaskpol1; }
    unsigned int antmaskpol2() { return m_hdr.m_antmaskpol2; }
    std::vector<std::string> antspol1() { return m_hdr.m_antspol1; }
    std::vector<std::string> antspol2() { return m_hdr.m_antspol2; }

    /** Beam tiling and steering parameters. **/
    int beamid() { return m_hdr.m_beamid; }
    int hostid() { return m_hdr.m_hostid; }
    int nbeams() { return m_hdr.m_nbeams; }
    std::string hostname() { return m_hdr.m_hostname; }
    int npcbaselines() { return m_hdr.m_npcbaselines; }
    int nbeamspernode() { return m_hdr.m_nbeamspernode; }
    std::vector<double> beamras() { return m_hdr.m_beamras; }
    std::vector<double> beamdecs() { return m_hdr.m_beamdecs; }

    /** Shared memory parameters **/
    int maxblks() { return FRBMAXBLKS; }
    int blksamps() { return FRBBLKSAMPS; }

    bool empty() { return m_bufptr->empty; }
    bool status() { return m_bufptr->status; }
    bool active() { return m_bufptr->active; }

    unsigned int curblk() { return m_bufptr->curblk; }
    unsigned int currec() { return m_bufptr->currec; }
    int begblk() { return (int)std::floor(curblk() / maxblks()) * maxblks(); }
    int endblk() { return begblk() + maxblks() - 1; }

    long blksize() { return blksamps() * nf(); }
    long size() { return maxblks() * blksize(); }

    double blktime() { return blksamps() * dt(); }
    double curtime() { return timeofblk(curblk()); }
    double begtime() { return timeofblk(begblk()); }
    double endtime() { return timeofblk(endblk()); }
    double timeofblk(int blk) { return blk * blktime(); }
    std::vector<std::chrono::system_clock::time_point> timestamps() { return m_timestamps; };

    /** Public methods. **/
    void open(MODE mode);
    void update();
    void close();

    void putblk(unsigned char* data, int beam, int blk);
    std::tuple<unsigned char*, size_t> getblk(int beam, int blk);
    void putblks(unsigned char* data, int beam, int blk0, int blkN);
    std::tuple<unsigned char*, size_t> getblks(int beam, int blk0, int blkN);
    std::tuple<unsigned char*, size_t> getslice(int beam, double tbeg, double tend);
    std::tuple<unsigned char*, size_t> getburst(int beam, double t0, double dm, double width);

  private:
    MODE m_mode;
    Header m_hdr;
    bool m_opened;

    /** Shared memory parameters. **/
    int m_hdrid;
    int m_bufid;
    BeamHeaderType* m_hdrptr;
    BeamBufferType* m_bufptr;
    unsigned char* m_dataptr;
    std::vector<std::chrono::system_clock::time_point> m_timestamps;

    /** Shared memory pointers. **/
    unsigned char* ptrtobeam(int beam);
    unsigned char* ptrtoblk(int beam, int blk);
    unsigned char* ptrtotime(int beam, double t);
  };
}  // namespace shazam

#endif
