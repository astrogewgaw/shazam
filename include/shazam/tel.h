#ifndef SHAZAM_TEL_H
#define SHAZAM_TEL_H

#include <sys/shm.h>

#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

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

  constexpr int TELHDRKEY = 1050;
  constexpr int TELBUFKEY = ShmKey;

  class TELRing {
  public:
    TELRing(MODE mode) {
      /** Open the header. **/
      m_mode = mode;

      switch (m_mode) {
        case READ: {
          /** Attach to header. **/
          m_hdrid = shmget(TELHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
          if (m_hdrid < 0) throw std::runtime_error("UNABLE TO GET HDR SHM ID. ABORT.");
          m_hdrptr = (BeamHeaderType*)shmat(m_hdrid, NULL, SHM_RDONLY);
          if ((void*)m_hdrptr == (void*)-1)
            throw std::runtime_error("FAILED TO LINK TO HDR SHM. ABORT");

          /** Read the header. **/
          ScanInfoType* scan = &(m_hdrptr->ScanTab[0]);

          /** Get some beam and host parameters early. **/
          m_beamid = m_hdrptr->BeamGenHdr.BeamHostID;
          m_hostid = m_hdrptr->BeamGenHdr.BeamHostID;
          m_hostname = m_hdrptr->BeamGenHdr.BeamHostName;

          /** Get data parameters. **/
          m_nbits = 8;
          m_nf = m_hdrptr->corr.corrpar.channels;
          m_fh = scan->source.freq[0] / 1e6;
          m_df = m_hdrptr->corr.corrpar.f_step / 1e6;
          m_flipped = scan->source.net_sign[0] == -1;
          m_dt = m_hdrptr->corr.daspar.gsb_final_bw * m_hdrptr->BeamGenHdr.SampInterval
                 / (m_hdrptr->corr.corrpar.clock);

          /** Some derived parameters. **/
          m_bw = m_nf * m_df;
          if (m_flipped) m_fh = m_fh + m_bw - 0.5 * m_df;
          m_fl = m_fh - m_bw + 0.5 * m_df;

          /** Get observation parameters. **/
          m_ra = scan->source.ra_app;
          m_dec = scan->source.dec_app;
          m_gtaccode = scan->proj.code;
          m_source = scan->source.object;
          m_gtactitle = scan->proj.title;
          m_observer = scan->proj.observer;
          m_nstokes = m_hdrptr->BeamGenHdr.NStokes[m_beamid];
          m_beammode = BEAMTYPES[m_hdrptr->BeamGenHdr.BeamType[m_beamid] - 1];

          /** Get antenna masks and antennas. **/
          unsigned int refantmask = 1;
          m_antmaskpol1 = m_hdrptr->BeamGenHdr.GAC_maskP1;
          for (int i = 0; i < 30; i++)
            if ((refantmask << i) & m_antmaskpol1) m_antspol1.push_back(ANTENNAS[i]);
          m_antmaskpol2 = m_hdrptr->BeamGenHdr.GAC_maskP2;
          for (int i = 0; i < 30; i++)
            if ((refantmask << i) & m_antmaskpol2) m_antspol2.push_back(ANTENNAS[i]);

          /** Get beam steering parameters. **/
          m_nbeams = m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeams;
          m_npcbaselines = m_hdrptr->BeamGenHdr.BeamSteeringParams.nPCBaselines;
          m_nbeamspernode = m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode;

          /** Get beam RA and DEC values. **/
          for (int i = 0; i < m_nbeamspernode; i++) {
            int b = m_beamid * m_nbeamspernode + i;
            m_beamras.push_back(m_hdrptr->BeamGenHdr.BeamSteeringParams.RA[b]);
            m_beamdecs.push_back(m_hdrptr->BeamGenHdr.BeamSteeringParams.DEC[b]);
          }

          int extrabuf = 64;
          int cursamps = 32 * 25;
          long int curtotalwords = cursamps * m_nf;
          long int currecsize = curtotalwords * WordSize / 2;

          long curshmdatasize = (MaxRecs + 1) * currecsize * m_nbeamspernode + extrabuf;
          curshmdatasize = curshmdatasize / PageSize + 1;
          curshmdatasize = curshmdatasize * PageSize;

          int shmdataoff = sizeof(GlobalInfoType) + extrabuf;
          shmdataoff = shmdataoff / PageSize + 1;
          shmdataoff = shmdataoff * PageSize;

          long curshmsize = curshmdatasize + shmdataoff;

          m_bufid = shmget(TELBUFKEY, curshmsize, SHM_RDONLY);
          if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET TEL SHM ID. ABORT.");
          m_bufptr = (GlobalInfoType*)shmat(m_bufid, NULL, SHM_RDONLY);
          if ((void*)m_bufptr == (void*)-1)
            throw std::runtime_error("FAILED TO OPEN TEL SHM. ABORT.");
          m_dataptr = (unsigned char*)m_bufptr;
          m_opened = true;
          break;
        }
        case WRITE: {
          /** Create (empty) header. **/
          m_hdrid = shmget(TELHDRKEY, sizeof(BeamHeaderType), IPC_CREAT | 0666);
          if (m_hdrid < 0) throw std::runtime_error("UNABLE TO GET HDR SHM ID. ABORT.");
          m_hdrptr = (BeamHeaderType*)shmat(m_hdrid, NULL, 0);
          if ((void*)m_hdrptr == (void*)-1)
            throw std::runtime_error("FAILED TO CREATE HDR SHM. ABORT");

          int extrabuf = 64;
          int cursamps = 32 * 25;
          long int curtotalwords = cursamps * m_nf;
          long int currecsize = curtotalwords * WordSize / 2;

          long curshmdatasize = (MaxRecs + 1) * currecsize * m_nbeamspernode + extrabuf;
          curshmdatasize = curshmdatasize / PageSize + 1;
          curshmdatasize = curshmdatasize * PageSize;

          int shmdataoff = sizeof(GlobalInfoType) + extrabuf;
          shmdataoff = shmdataoff / PageSize + 1;
          shmdataoff = shmdataoff * PageSize;

          long curshmsize = curshmdatasize + shmdataoff;

          m_bufid = shmget(TELBUFKEY, curshmsize, IPC_CREAT | 0777);
          if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET TEL SHM ID. ABORT.");
          m_bufptr = (GlobalInfoType*)shmat(m_bufid, NULL, 0);
          if ((void*)m_bufptr == (void*)-1)
            throw std::runtime_error("FAILED TO OPEN TEL SHM. ABORT.");
          m_dataptr = (unsigned char*)m_bufptr;

          /** Get timestamps. **/
          for (int ii = 0; ii < maxblks(); ++ii) {
            m_timestamps.push_back(std::chrono::system_clock::time_point{
                std::chrono::seconds{m_bufptr->rec[ii].timestamp_gps.tv_sec}
                + std::chrono::microseconds{m_bufptr->rec[ii].timestamp_gps.tv_usec}
                + std::chrono::nanoseconds{(long)m_bufptr->rec[ii].blk_nano}});
          }

          m_opened = true;
          break;
        }
      }
    }

    ~TELRing() {}

    MODE mode() { return m_mode; }
    bool opened() { return m_opened; }

    /** Data parameters. **/
    int nf() { return m_nf; }
    double fh() { return m_fh; }
    double fl() { return m_fl; }
    double df() { return m_df; }
    double bw() { return m_bw; }
    double dt() { return m_dt; }
    int nbits() { return m_nbits; }
    int nstokes() { return m_nstokes; }
    bool flipped() { return m_flipped; }

    /** Observation parameters. **/
    double ra() { return m_ra; }
    double dec() { return m_dec; }
    std::string source() { return m_source; }
    std::string beammode() { return m_beammode; }
    std::string observer() { return m_observer; }
    std::string gtaccode() { return m_gtaccode; }
    std::string gtactitle() { return m_gtactitle; }
    unsigned int antmaskpol1() { return m_antmaskpol1; }
    unsigned int antmaskpol2() { return m_antmaskpol2; }
    std::vector<std::string> antspol1() { return m_antspol1; }
    std::vector<std::string> antspol2() { return m_antspol2; }

    /** Beam tiling and steering parameters. **/
    int beamid() { return m_beamid; }
    int hostid() { return m_hostid; }
    int nbeams() { return m_nbeams; }
    std::string hostname() { return m_hostname; }
    int npcbaselines() { return m_npcbaselines; }
    int nbeamspernode() { return m_nbeamspernode; }
    std::vector<double> beamras() { return m_beamras; }
    std::vector<double> beamdecs() { return m_beamdecs; }

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
    bool m_opened;

    /** Data parameters. **/
    int m_nf;
    double m_fh;
    double m_fl;
    double m_df;
    double m_bw;
    double m_dt;
    int m_nbits;
    int m_nstokes;
    bool m_flipped;

    /** Observation parameters. **/
    double m_ra;
    double m_dec;
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
    std::string m_hostname;
    int m_npcbaselines;
    int m_nbeamspernode;
    std::vector<double> m_beamras;
    std::vector<double> m_beamdecs;

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
