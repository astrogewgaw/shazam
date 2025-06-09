#ifndef MULTIHDR_H
#define MULTIHDR_H

#ifdef __cplusplus
extern "C" {
#endif

enum { USB_130, USB_175, LSB_130, LSB_175, MAX_BANDS };
enum { RRLL, RRRL, RR__, MACMODES };
enum { DAS_CARD0, DAS_CARD1, DAS_CARDS, BEDAS_CARDS, GPUDAS_CARDS };
enum {
  IndianPolar,
  UsbPolar,
  LsbPolar,
  UsbCopy,
  LsbCopy,
  AllU130,
  AllU175,
  AllL130,
  AllL175,
  arar_arar,
  alal_alal,
  brbr_brbr,
  blbl_blbl,
  aral_brbl,
  aral_alar,
  brbl_blbr,
  arbr_albl,
  DpcMuxVals
};

enum { NAMELEN = 32, DATELEN = 32 };

enum {
  MAX_ANTS = 32,
  MAX_SAMPS = 64,
  MAX_FFTS = MAX_SAMPS,
  MAC_CARDS = 33,
  MAX_BASE = DAS_CARDS * 32 * MAC_CARDS,
  MAX_CHANS = 32768,
  POLS = 2
};

enum { MAX_PROJECTS = 500, MAX_SCANS = 512, MAX_BM_STEER = 2568 };
enum { MAX_ARRAYS = 1 };
enum { MAX_BEAMS = 4 };
enum { MAX_GATES = 16 };
enum { LittleEndian = 1, BigEndian = 0 };
enum { TransitObs = 32768 };

enum {
  TimeSize = sizeof(double),
  WtSize = 2 * sizeof(float),
  ActiveScans = 64,
  DataFlagSize = ActiveScans * sizeof(int)
};

typedef struct interfaces {
  char iface[30];
} Interface;

typedef struct {
  char name[4];
  unsigned char samp_id[MAX_BANDS];
  double bx;
  double by;
  double bz;
  double d0[MAX_BANDS];
  double p0[MAX_BANDS];
} AntennaParType;

typedef struct {
  int macs;
  int channels;
  int pols;
  int sta;
  int cntrl;
  int statime;
  int iabeam;
  int pabeam1;
  int pabeam2;
  float iabeam_res;
  float pabeam1_res;
  float pabeam2_res;
  float f_step;
  float clock;
  unsigned char dpcmux;
  unsigned char clksel;
  unsigned char fftmode;
  unsigned char macmode;
} CorrParType;

typedef struct {
  int antmask;
  int samplers;
  int baselines;
  int channels;
  int lta;
  int gsb_maxchan;
  int gsb_fstop;
  short bandmask;
  short mode;
  short gsb_stokes;
  unsigned short chan_num[MAX_CHANS];
  double mjd_ref;
  double t_unit;
  double gsb_acq_bw;
  double gsb_final_bw;
  double dut1;
} DasParType;

typedef struct {
  char object[NAMELEN];
  struct {
    float i;
    float q;
    float u;
    float v;
  } flux;
  double mjd0;
  double ra_app;
  double dec_app;
  double ra_mean;
  double dec_mean;
  double dra;
  double ddec;
  double freq[2];
  double first_lo[2];
  double bb_lo[2];
  double rest_freq[2];
  double lsrvel[2];
  double ch_width;
  int id;
  int net_sign[MAX_BANDS];
  int mode;
  int dum1;
  unsigned int antmask;
  unsigned short bandmask;
  unsigned short dum2;
  short calcode;
  short qual;
} SourceParType;

typedef struct {
  char object[NAMELEN];
  char date_obs[DATELEN];
  char observer[NAMELEN];
  char project[NAMELEN];
  char code[8];
  double ra_date;
  double dec_date;
  double mjd_ref;
  double dra;
  double ddec;
  double period;
  int antmask;
  int bandmask;
  int flag;
  int seq;
  int fold_flg;
  float integ;
  float f_step;
  float dm;
  int rf[2];
  int first_lo[2];
  int bb_lo[2];
  float BW;
  int ref_ch;
  int i_side_band;
} ScanParType;

typedef struct {
  char code[8];
  char observer[NAMELEN];
  char title[NAMELEN];
  unsigned int antmask;
  unsigned short bandmask;
  unsigned short seq;
} ProjectType;

typedef struct {
  int status;
  float t;
  ProjectType proj;
  SourceParType source;
} ScanInfoType;

typedef struct {
  unsigned char ant_id;
  unsigned char band;
  unsigned char fft_id;
  unsigned char dpc;
} SamplerType;

typedef struct {
  SamplerType samp[2];
  unsigned char card;
  unsigned char chip;
  unsigned char pol;
  unsigned char word_incr;
} BaseParType;

typedef struct {
  unsigned char endian;
  unsigned char dummy[7];
  char version[NAMELEN];
  char bandname[MAX_BANDS][8];
  AntennaParType antenna[MAX_ANTS];
  SamplerType sampler[MAX_SAMPS];
  BaseParType baseline[MAX_BASE];
  CorrParType corrpar;
  DasParType daspar;
} CorrType;

typedef struct {
  unsigned int ext_model;
  unsigned int idle;
  unsigned int stop;
  unsigned int userflag;
} AntennaFlagType;

typedef struct {
  float phase;
  float dp;
  float delay;
  float dd;
} ModelParType;

typedef struct {
  float phase;
  float dp;
  float delay;
  float dslope;
} DataParType;

typedef struct {
  double t0;
  int ant_id;
  int band;
  ModelParType par;
} ExtModelType;

typedef struct {
  int BitPackFlag;
  int DataBitFormat;
  int DataLSBBits;
  int DataBitWinSize;
} BitPackType;

typedef struct {
  int FoldingFlag;
  double FoldingPeriod;
} FoldingType;

typedef struct {
  int BBRFIFilterFlag;
  char ReplacementType[64];
  int ReplaceValue;
  float Threshold;
  char Strategy[32];
  char param1[64];
  char param2[64];
  char param3[64];
} BBRFIFilterType;

typedef struct {
  int BandEqFlag;
  float optlevel;
  int beg_chan;
  int end_chan;
  float rollOff;
  int taps;
  char param1[64];
  char param2[64];
  char param3[64];
} BandEqType;

typedef struct {
  int BeamSteeringFlag;
  int nPCBaselines;
  int nBeamHosts;
  int nSteeringBeams;
  short int nSteeringBeamsPerNode;
  float RA[MAX_BM_STEER];
  float DEC[MAX_BM_STEER];
  int Beam_index[MAX_BM_STEER];
  int Beam_subindex[MAX_BM_STEER];
  char param2[64];
  char param3[64];
  char param4[64];
} BeamSteeringType;

typedef struct {
  char Feature1Title[64];
  int Feature1Flag;
  char param1[64];
  char param2[64];
  char param3[64];
  char param4[64];
} Feature1Type;

typedef struct {
  char Feature2Title[64];
  int Feature2Flag;
  char param1[64];
  char param2[64];
  char param3[64];
  char param4[64];
} Feature2Type;

typedef struct {
  char Feature3Title[64];
  int Feature3Flag;
  char param1[64];
  char param2[64];
  char param3[64];
  char param4[64];
} Feature3Type;

typedef struct {
  char Feature4Title[64];
  int Feature4Flag;
  char param1[64];
  char param2[64];
  char param3[64];
  char param4[64];
} Feature4Type;

typedef struct {
  int WalshFlag;
  char param1[64];
  char param2[64];
  char param3[64];
} WalshType;

typedef struct {
  int PFBFlag;
  int Taps;
  char param1[64];
  char param2[64];
  char param3[64];
} PFBType;

typedef struct {
  int CDFlag;
  unsigned int in_channels;
  unsigned int out_subbands;
  float freq;
  float bw;
  float DM;
  int sb_flag;
  int vlt_bits;
  int write_int;
  float sampling_resolution;
} CDType;

typedef struct {
  int DedispersionFlag;
  double dD_RF;
  double dD_DM;
  double dD_BW;
  int dD_RefChan;
  int dD_iSB;
} dDType;

typedef struct {
  int GatingFlag;
  int NoOfGates;
  int GtPeriod;
  int GtNoOfBins;
  int GtStartBinNo[MAX_GATES];
  int GtStopBinNo[MAX_GATES];
} GatingType;

typedef struct {
  int FilterFlag;
  int FftLength;
  int NBands;
  int FilterDataType;
} RfiFilterType;

typedef struct {
  int BeamHostID;
  char BeamHostName[128];
  int PreTimeInt;
  double SampInterval;
  int GAC_maskP1;
  int GAC_maskP2;
  unsigned long GAC_maskFp[MAX_BEAMS];
  int GAC_vlt_maskP1[MAX_BEAMS];
  int GAC_vlt_maskP2[MAX_BEAMS];
  int PostFreqInt[MAX_BEAMS];
  int PostTimeInt[MAX_BEAMS];
  double PrePeriod[MAX_BEAMS];
  float DM[MAX_BEAMS];
  int BeamType[MAX_BEAMS];
  int NStokes[MAX_BEAMS];
  int OutputDataFormat;
  BitPackType BitPackParams;
  dDType dDParams;
  CDType CDParams;
  FoldingType FoldingParams;
  GatingType GatingParams;
  RfiFilterType RfiFilterParams;
  PFBType PFBparams;
  WalshType WalshParams;
  BandEqType BandEqParams;
  BBRFIFilterType BBRFIFilterParams;
  BeamSteeringType BeamSteeringParams;
  Feature1Type Feature1Params;
  Feature2Type Feature2Params;
  Feature3Type Feature3Params;
  Feature4Type Feature4Params;
} BeamGenType;

typedef struct {
  short int SyncCmd;
  short int VersionMajor;
  short int VersionMinor;
  ProjectType Project[MAX_PROJECTS];
  ScanInfoType ScanTab[MAX_SCANS];
  ScanInfoType ScanTab1[MAX_SCANS];
  ScanInfoType ScanTab2[MAX_SCANS];
  ScanInfoType ScanTab3[MAX_SCANS];
  CorrType corr;
  ScanParType scanpar[MAX_ARRAYS];
  BeamGenType BeamGenHdr;
} BeamHeaderType;

#define AntennaTypeSize sizeof(AntennaType)
#define MacFftTypeSize sizeof(MacFftType)
#define SamplerTypeSize sizeof(SamplerType)
#define DataParTypeSize sizeof(DataParType)
#define CorrTypeSize sizeof(CorrType)
#define CorrSize sizeof(CorrType)
#define Corr2Size 16192

typedef struct {
  int in0;
  int in1;
  int out0;
  int out1;
} IndexType;

typedef struct {
  double bxcd;
  double bycd;
  double bzsd;
  double fixed;
  double phase;
  double bb_lo;
  double freq;
  double f_step;
} GeomType;

#define DAS_H_KEY 1030
#define DAS_D_KEY 1031
#define DAS_H0_KEY 1032
#define DAS_D0_KEY 1033
#define DAS_HDRSIZE 200000
#define DAS_BUFSIZE 545259520

typedef struct {
  int s0;
  int s1;
  int card;
  int delay;
  int p0;
  int pstep0;
  int p1;
  int pstep1;
  int fstc;
  float p2fstc;
  float fstc_step;
} FftDelayParType;

typedef struct {
  double clock;
  double t_update;
  double pc_time;
  double t_off;
  double delay_step;
  double fstc_scale;
  double nco_tick_time;
  int cycle;
  int seq;
  int cards;
  unsigned char dpcmux;
  unsigned char clksel;
  unsigned char fftmode;
  unsigned char macmode;
  ModelParType par[MAX_SAMPS];
  FftDelayParType fdpar[MAX_SAMPS / 2];
} ModelInfoType;

typedef struct {
  int active;
  int status;
  int scan;
  int scan_off;
  CorrType corr;
  ModelInfoType model;
  BeamHeaderType BeamHeader;
  char buf[DAS_HDRSIZE];
} DasHdrType;

enum {
  BufMarked = 1,
  BufReady = 1 << 1,
  Rack0Bad = 1 << 2,
  Rack1Bad = 1 << 3,
  Rack2Bad = 1 << 4,
  Rack3Bad = 1 << 5,
  Rack4Bad = 1 << 6,
  Rack5Bad = 1 << 7,
  BufFinish = 1 << 8,
  MaxDataBuf = 100
};

enum { MAX_EVENTS = 50000 };

typedef struct {
  float t;
  unsigned char type;
  unsigned char cmd;
  unsigned short seq;
  int flag_num;
  int scan_num;
} EventLogType;

typedef struct {
  int t_off;
  int wt_off;
  int par_off;
  int data_off;
  int data_words;
  int par_words;
  int wordsize;
} RecOffsetType;

typedef struct {
  int off;
  BaseParType base;
  char name[12];
} MacWordType;

typedef struct {
  unsigned char seq;
  unsigned char scan;
  int status;
  int recl;
  int seqnum;
  RecOffsetType off;
  MacWordType *mac;
  char *buf;
} RecBufType;

typedef struct {
  int active;
  int status;
  unsigned short events;
  unsigned short flags;
  unsigned short starts;
  unsigned short stops;
  CorrType corr;
  AntennaFlagType flag[MAX_EVENTS][MAX_BANDS];
  EventLogType event[MAX_EVENTS];
  ScanInfoType scaninfo[MAX_SCANS];
  RecOffsetType offset;
} DataInfoType;

typedef struct {
  int flag;
  int rec;
  int seqnum;
  unsigned short flag_seq;
  unsigned short newstate;
} DataTabType;

typedef struct {
  int flag;
  int blocksize;
  int maxblocks;
  int cur_block;
  int first_block;
  int cur_rec;
  DataTabType dtab[MaxDataBuf];
  char buf[DAS_BUFSIZE];
} DataBufType;

#define CHANNEL 4096
#define FFTLEN (2 * CHANNEL)
#define NCHAN 64
#define CORRLEN (NCHAN * FFTLEN)
#define NUM_ANT 1
#define NUM_ACQ 1

#define BandWidth 400000000
#define BLOCKTIME (8 / 32)
#define SAMP_CLK 0.000000005

#ifdef __cplusplus
}
#endif

#include <stdexcept>
#include <string>
#include <sys/shm.h>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

constexpr int MULTIHDRKEY = 1050;

constexpr char ANTENNAS[30][4] = {
    "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C08", "C09", "C10",
    "C11", "C12", "C13", "C14", "E02", "E03", "E04", "E05", "E06", "S01",
    "S02", "S03", "S04", "S06", "W01", "W02", "W03", "W04", "W05", "W06"};

constexpr char BEAMTYPES[7][6] = {"IA",  "PA",   "VLT", "PC",
                                  "CDP", "PASV", "MISC"};

class MultiHeader {
public:
  MultiHeader() { link(); };
  ~MultiHeader() { unlink(); };

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

  void link();
  void unlink();
  nb::dict asdict();

  friend class MultiTELSHM;
  friend class MultiFRBSHM;

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
  bool linked;
  BeamHeaderType *m_hdrptr;
};

void initmultihdr(nb::module_ m);

#endif
