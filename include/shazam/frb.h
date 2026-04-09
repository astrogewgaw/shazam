#ifndef SHAZAM_FRB_H
#define SHAZAM_FRB_H

#include <sys/shm.h>

#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "hdr.h"

namespace shazam {
  constexpr int FRBMAXBLKS = 12;
  constexpr int FRBHDRKEY = 2031;
  constexpr int FRBBUFKEY = 2032;
  constexpr int FRBFFTSAMPS = 800;
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
    FRBRing(MODE mode) {
      /** Open the header. **/
      m_mode = mode;
      m_opened = false;

      switch (m_mode) {
        case READ: {
          /** Attach to header. **/
          m_hdrid = shmget(FRBHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
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

          /** Calculate size of buffer. **/
          long BLKSIZE = (long)FRBBLKSAMPS * (long)(m_nf);
          long BUFSIZE = BLKSIZE * (long)FRBMAXBLKS * (long)(m_nbeamspernode);
          long FRBSHMSIZE = sizeof(BeamBufferType) + BUFSIZE;

          /** Attach to the buffer. **/
          m_bufid = shmget(FRBBUFKEY, FRBSHMSIZE, SHM_RDONLY);
          if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET FRB SHM ID. ABORT.");
          m_bufptr = (BeamBufferType*)shmat(m_bufid, NULL, SHM_RDONLY);
          if ((void*)m_bufptr == (void*)-1)
            throw std::runtime_error("FAILED TO OPEN FRB SHM. ABORT.");
          m_dataptr = ((unsigned char*)m_bufptr) + sizeof(BeamBufferType);

          /** Get timestamps. **/
          for (int ii = 0; ii < maxblks(); ++ii) {
            m_timestamps.push_back(std::chrono::system_clock::time_point{
                std::chrono::seconds{m_bufptr->timestamps[ii].tv_sec}
                + std::chrono::microseconds{m_bufptr->timestamps[ii].tv_usec}
                + std::chrono::nanoseconds{(long)m_bufptr->nanoseconds[ii]}});
          }

          m_opened = true;
          break;
        }
        case WRITE: {
          /** Create (empty) header. **/
          m_hdrid = shmget(FRBHDRKEY, sizeof(BeamHeaderType), IPC_CREAT | 0666);
          if (m_hdrid < 0) throw std::runtime_error("UNABLE TO GET HDR SHM ID. ABORT.");
          m_hdrptr = (BeamHeaderType*)shmat(m_hdrid, NULL, 0);
          if ((void*)m_hdrptr == (void*)-1)
            throw std::runtime_error("FAILED TO CREATE HDR SHM. ABORT");

          ScanInfoType* scan = &(m_hdrptr->ScanTab[0]);

          /** Set some beam and host parameters early. **/
          m_hdrptr->BeamGenHdr.BeamHostID = m_beamid;
          m_hdrptr->BeamGenHdr.BeamHostID = m_hostid;
          strcpy(m_hdrptr->BeamGenHdr.BeamHostName, m_hostname.c_str());

          /** Set data parameters. **/
          scan->source.freq[0] = m_fh * 1e6;
          m_hdrptr->corr.corrpar.channels = m_flipped ? m_fl : m_nf;
          m_hdrptr->corr.daspar.gsb_final_bw = 1;
          m_hdrptr->BeamGenHdr.PostTimeInt[0] = 1;
          m_hdrptr->BeamGenHdr.PostFreqInt[0] = 1;
          m_hdrptr->corr.daspar.gsb_acq_bw = m_bw;
          m_hdrptr->corr.corrpar.f_step = m_df * 1e6;
          scan->source.net_sign[0] = m_flipped ? -1 : 1;
          m_hdrptr->corr.corrpar.clock = 2.0 * m_bw * 1e6;
          m_hdrptr->BeamGenHdr.SampInterval
              = m_dt * m_hdrptr->corr.corrpar.clock / m_hdrptr->corr.daspar.gsb_final_bw;

          /** Set observation parameters. **/
          scan->source.ra_app = m_ra;
          scan->source.dec_app = m_dec;
          strcpy(scan->proj.code, m_gtaccode.c_str());
          strcpy(scan->source.object, m_source.c_str());
          strcpy(scan->proj.title, m_gtactitle.c_str());
          strcpy(scan->proj.observer, m_observer.c_str());
          m_hdrptr->BeamGenHdr.NStokes[m_beamid] = m_nstokes;
          if (m_beammode == "IA") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 0;
          } else if (m_beammode == "PA") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 1;
          } else if (m_beammode == "VLT") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 2;
          } else if (m_beammode == "PC") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 3;
          } else if (m_beammode == "CDP") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 4;
          } else if (m_beammode == "PASV") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 5;
          } else if (m_beammode == "MISC") {
            m_hdrptr->BeamGenHdr.BeamType[m_beamid] = 6;
          }

          /** Set antenna masks and antennas. **/
          m_hdrptr->BeamGenHdr.GAC_maskP1 = m_antmaskpol1;
          m_hdrptr->BeamGenHdr.GAC_maskP2 = m_antmaskpol2;

          /** Set beam steering parameters. **/
          m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeams = m_nbeams;
          m_hdrptr->BeamGenHdr.BeamSteeringParams.nPCBaselines = m_npcbaselines;
          m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode = m_nbeamspernode;

          /** Get beam RA and DEC values. **/
          for (int i = 0; i < m_nbeamspernode; i++) {
            int b = m_beamid * m_nbeamspernode + i;
            m_hdrptr->BeamGenHdr.BeamSteeringParams.RA[b] = m_beamras[i];
            m_hdrptr->BeamGenHdr.BeamSteeringParams.DEC[b] = m_beamdecs[i];
          }

          /** Calculate size of buffer. **/
          long BLKSIZE = (long)FRBBLKSAMPS * (long)(m_nf);
          long BUFSIZE = BLKSIZE * (long)FRBMAXBLKS * (long)(m_nbeamspernode);
          long FRBSHMSIZE = sizeof(BeamBufferType) + BUFSIZE;

          m_bufid = shmget(FRBBUFKEY, FRBSHMSIZE, IPC_CREAT | 0666);
          if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET FRB SHM ID. ABORT.");
          m_bufptr = (BeamBufferType*)shmat(m_bufid, NULL, 0);
          if ((void*)m_bufptr == (void*)-1)
            throw std::runtime_error("FAILED TO OPEN FRB SHM. ABORT.");
          m_dataptr = ((unsigned char*)m_bufptr) + sizeof(BeamBufferType);

          m_opened = true;
          break;
        }
      }
    }

    ~FRBRing() {}

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
