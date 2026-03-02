#ifndef SHAZAM_TEL_H
#define SHAZAM_TEL_H

#include <sys/shm.h>

#include <chrono>
#include <cmath>
#include <string>
#include <tuple>

#include "hdr.h"

namespace shazam {
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

  enum { MaxGPS = 16, MaxBLK = 64, MaxRecs = 8, ExtraWords = 32, ExtraBuf = ExtraWords * WordSize };

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
    unsigned short* begp;
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
    void* shmp;
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
    ShmDataSize = (((MaxRecs + 1) * RecSize + ExtraBuf) / PageSize + 1) * PageSize,
    ShmSize = ShmDataOff + ShmDataSize
  };

#define FFT_CYCLE CHANNEL * 2
#define BLK_INTEG 8192
#define BASE_CLK 400.0e6

#ifdef __cplusplus
  }
#endif

  class TELRing {
  public:
    TELRing()
        : m_hdr(),
          m_hdrid(0),
          m_bufid(0),
          m_mode(READ),
          m_hdrptr(NULL),
          m_bufptr(NULL),
          m_opened(false),
          m_dataptr(NULL) {}

    ~TELRing() {}

    MODE mode() { return m_mode; }
    Header header() { return m_hdr; }
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
    int maxblks() { return MaxRecs; }
    int blksamps() { return 32 * 25; }

    bool acqover() { return m_bufptr->acq_flag & AcqOver; }
    bool gpsok() { return m_bufptr->acq_flag & GPSpresent; }
    bool acqok() { return !(m_bufptr->acq_flag & UnInitialized); }

    bool marked() { return m_bufptr->rec[currec()].rec_flag & Marked; }
    bool dataok() { return m_bufptr->rec[currec()].rec_flag & GoodData; }
    bool blkok() { return !(m_bufptr->rec[currec()].rec_flag & BlockErr); }
    bool timeok() { return !(m_bufptr->rec[currec()].rec_flag & TimeErr); }
    bool noinit() { return m_bufptr->rec[currec()].rec_flag & UnInitialized; }

    unsigned int curblk() { return m_bufptr->rec[currec()].rec_seq; }
    unsigned int currec() { return (m_bufptr->rec_ind - 1) % maxblks(); }
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

  private:
    MODE m_mode;
    Header m_hdr;
    bool m_opened;

    /** Shared memory parameters. **/
    int m_hdrid;
    int m_bufid;
    BeamHeaderType* m_hdrptr;
    GlobalInfoType* m_bufptr;
    unsigned char* m_dataptr;
    std::vector<std::chrono::system_clock::time_point> m_timestamps;

    /** Shared memory pointers. **/
    unsigned char* ptrtobeam(int beam);
    unsigned char* ptrtoblk(int beam, int blk);
    unsigned char* ptrtotime(int beam, double t);
  };
}  // namespace shazam

#endif
