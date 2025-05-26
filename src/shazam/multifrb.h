#ifndef MULTIFRB_H
#define MULTIFRB_H

#include <cstring>
#include <ctime>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/shm.h>
#include <sys/time.h>
#include <vector>

#include "multihdr.h"

constexpr int FRBHDRKEY = 2031;
constexpr int FRBBUFKEY = 2032;

constexpr int FFTSAMPS = 800;
constexpr int FRBMAXBLKS = 12;
constexpr int FRBBLKSAMPS = FFTSAMPS * 32;

constexpr char ANTLIST[30][4] = {
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
std::string safe_strfmt(const std::string &format, Args... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;
  if (size_s <= 0) throw std::runtime_error("Error during formatting.");
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
  MultiFRBSHM() {
    /** Attach to header. **/
    int HDRID = shmget(FRBHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
    if (HDRID < 0) throw std::runtime_error("");
    hdrptr = (BeamHeaderType *)shmat(HDRID, 0, SHM_RDONLY);
    if ((void *)hdrptr == (void *)-1) throw std::runtime_error("");

    /** Read in all header parameters... **/
    ScanInfoType *scan = &(hdrptr->ScanTab[0]);

    /** Get some beam and host parameters early. **/
    m_beamid = hdrptr->BeamGenHdr.BeamHostID;
    m_hostid = hdrptr->BeamGenHdr.BeamHostID;
    m_hostname = hdrptr->BeamGenHdr.BeamHostName;

    /** Get data parameters. **/
    m_nbits = 8;
    m_nf = hdrptr->corr.corrpar.channels;
    m_fh = scan->source.freq[0] / 1e6;
    m_df = hdrptr->corr.corrpar.f_step / 1e6;
    m_flipped = scan->source.net_sign[0] == -1;
    m_dt = hdrptr->corr.daspar.gsb_final_bw * hdrptr->BeamGenHdr.SampInterval /
           (hdrptr->corr.corrpar.clock);

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
    m_nstokes = hdrptr->BeamGenHdr.NStokes[m_beamid];
    m_beammode = BEAMTYPES[hdrptr->BeamGenHdr.BeamType[m_beamid] - 1];

    /** Get antennas mask and antennas. **/
    unsigned int refantmask = 1;
    m_antmaskpol1 = hdrptr->BeamGenHdr.GAC_maskP1;
    for (int i = 0; i < 30; i++)
      if ((refantmask << i) & m_antmaskpol1) m_antspol1.push_back(ANTLIST[i]);
    m_antmaskpol2 = hdrptr->BeamGenHdr.GAC_maskP2;
    for (int i = 0; i < 30; i++)
      if ((refantmask << i) & m_antmaskpol2) m_antspol2.push_back(ANTLIST[i]);

    /** Get beam steering parameters. **/
    m_nbeams = hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeams;
    m_npcbaselines = hdrptr->BeamGenHdr.BeamSteeringParams.nPCBaselines;
    m_nbeamspernode =
        hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode;

    /** Get beam RA and DEC values. **/
    for (int i = 0; i < m_nbeamspernode; i++) {
      int b = m_beamid * m_nbeamspernode + i;
      m_beamras.push_back(hdrptr->BeamGenHdr.BeamSteeringParams.RA[b]);
      m_beamdecs.push_back(hdrptr->BeamGenHdr.BeamSteeringParams.DEC[b]);
    }

    /** Calculate size of buffer. **/
    long BLKSIZE = (long)FRBBLKSAMPS * (long)m_nf;
    long BUFSIZE = BLKSIZE * (long)FRBMAXBLKS * (long)m_nbeamspernode;
    long FRBSHMSIZE = sizeof(BeamBufferType) + BUFSIZE;

    /** Attach to buffer. **/
    int BUFID = shmget(FRBBUFKEY, FRBSHMSIZE, SHM_RDONLY);
    if (BUFID < 0) throw std::runtime_error("");
    bufptr = (BeamBufferType *)shmat(BUFID, 0, SHM_RDONLY);
    if ((void *)bufptr == (void *)-1) throw std::runtime_error("");
    dataptr = ((unsigned char *)bufptr) + sizeof(BeamBufferType);

    /** Get index of current record and block being written. **/
    if (bufptr->empty) {
      curblk = bufptr->curblk;
      currec = bufptr->currec;
    } else {
      curblk = bufptr->curblk;
      currec = (bufptr->currec - 1) % FRBMAXBLKS;
    }

    /** Get date and time of observation. **/
    double seconds;
    struct tm *localt;
    double nanoseconds = 0;
    struct timeval timestamp;

    memcpy(&timestamp, &(bufptr->timestamps[currec]), sizeof(struct timeval));
    memcpy(&nanoseconds, &(bufptr->nanoseconds[currec]), sizeof(nanoseconds));

    localt = localtime(&timestamp.tv_sec);
    seconds = localt->tm_sec;
    seconds += (timestamp.tv_usec) / 1e6;
    seconds += (nanoseconds) / 1e9;

    m_obsdate = safe_strftime("%d/%m/%Y", localt);
    m_obstime = safe_strfmt("%02d:%02d:%012.9lf", localt->tm_hour,
                            localt->tm_min, seconds);
  };

  ~MultiFRBSHM() {};

  int nf() { return m_nf; };
  double fh() { return m_fh; };
  double fl() { return m_fl; };
  double df() { return m_df; };
  double bw() { return m_bw; };
  double dt() { return m_dt; };
  int nbits() { return m_nbits; };
  int nstokes() { return m_nstokes; };
  bool flipped() { return m_flipped; };

  double ra() { return m_ra; };
  double dec() { return m_dec; };
  int beamid() { return m_beamid; };
  int hostid() { return m_hostid; };
  int nbeams() { return m_nbeams; };
  std::string source() { return m_source; };
  std::string obsdate() { return m_obsdate; };
  std::string obstime() { return m_obstime; };
  int npcbaselines() { return m_npcbaselines; };
  int nbeamspernode() { return m_nbeamspernode; };
  std::string hostname() { return m_hostname; };
  std::string beammode() { return m_beammode; };
  std::string observer() { return m_observer; };
  std::string gtaccode() { return m_gtaccode; };
  std::string gtactitle() { return m_gtactitle; };
  unsigned int antmaskpol1() { return m_antmaskpol1; };
  unsigned int antmaskpol2() { return m_antmaskpol2; };
  std::vector<double> beamras() { return m_beamras; };
  std::vector<double> beamdecs() { return m_beamdecs; };
  std::vector<std::string> antspol1() { return m_antspol1; };
  std::vector<std::string> antspol2() { return m_antspol2; };

private:
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

  double m_ra;
  double m_dec;
  int m_beamid;
  int m_hostid;
  int m_nbeams;
  int m_npcbaselines;
  int m_nbeamspernode;
  std::string m_obsdate;
  std::string m_obstime;
  std::string m_source;
  std::string m_hostname;
  std::string m_beammode;
  std::string m_observer;
  std::string m_gtaccode;
  std::string m_gtactitle;
  unsigned int m_antmaskpol1;
  unsigned int m_antmaskpol2;
  std::vector<double> m_beamras;
  std::vector<double> m_beamdecs;
  std::vector<std::string> m_antspol1;
  std::vector<std::string> m_antspol2;

  int currec;
  int curblk;
  BeamHeaderType *hdrptr;
  BeamBufferType *bufptr;
  unsigned char *dataptr;
};

#endif
