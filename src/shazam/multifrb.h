#ifndef MULTIFRB
#define MULTIFRB

#include <nanobind/nanobind.h>
#include <sys/shm.h>
#include <sys/time.h>
#include <unistd.h>

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
  double bx, by, bz; /* metres */
  double d0[MAX_BANDS], p0[MAX_BANDS];
} AntennaParType;

typedef struct {
  int macs, channels, pols, sta, cntrl, statime, iabeam, pabeam1, pabeam2;
  float iabeam_res, pabeam1_res, pabeam2_res, f_step, clock;
  unsigned char dpcmux, clksel, fftmode, macmode;
} CorrParType;

typedef struct {
  int antmask, samplers, baselines, channels, lta, gsb_maxchan, gsb_fstop;
  short bandmask, mode, gsb_stokes;
  unsigned short chan_num[MAX_CHANS];
  double mjd_ref, t_unit, gsb_acq_bw, gsb_final_bw;
  double dut1;
} DasParType;

typedef struct {
  char object[NAMELEN];
  struct {
    float i, q, u, v;
  } flux;
  double mjd0;
  double ra_app, dec_app, ra_mean, dec_mean, dra, ddec;
  double freq[2], first_lo[2], bb_lo[2];
  double rest_freq[2], lsrvel[2];
  double ch_width;
  int id, net_sign[MAX_BANDS], mode, dum1;
  unsigned int antmask;
  unsigned short bandmask, dum2;
  short calcode, qual;
} SourceParType;

typedef struct {
  char object[NAMELEN], date_obs[DATELEN];
  char observer[NAMELEN], project[NAMELEN], code[8];
  double ra_date, dec_date, mjd_ref, dra, ddec, period;
  int antmask, bandmask, flag, seq, fold_flg;
  float integ, f_step, dm;
  int rf[2], first_lo[2], bb_lo[2];
  float BW;
  int ref_ch, i_side_band;
} ScanParType;

typedef struct {
  char code[8], observer[NAMELEN], title[NAMELEN];
  unsigned int antmask;
  unsigned short bandmask, seq;
} ProjectType;

typedef struct {
  int status;
  float t;
  ProjectType proj;
  SourceParType source;
} ScanInfoType;

typedef struct {
  unsigned char ant_id, band, fft_id, dpc;
} SamplerType;

typedef struct {
  SamplerType samp[2];
  unsigned char card, chip, pol, word_incr;
} BaseParType;

typedef struct {
  unsigned char endian, dummy[7];
  char version[NAMELEN];
  char bandname[MAX_BANDS][8];
  AntennaParType antenna[MAX_ANTS];
  SamplerType sampler[MAX_SAMPS];
  BaseParType baseline[MAX_BASE];
  CorrParType corrpar;
  DasParType daspar;
} CorrType;

typedef struct {
  unsigned int ext_model, idle, stop, userflag;
} AntennaFlagType;

typedef struct {
  float phase, dp, delay, dd;
} ModelParType;
typedef struct {
  float phase, dp, delay, dslope;
} DataParType;

typedef struct {
  double t0;
  int ant_id, band;
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
  int in0, in1, out0, out1;
} IndexType;

typedef struct {
  double bxcd, bycd, bzsd, fixed, phase, bb_lo, freq, f_step;
} GeomType;

#define DAS_H_KEY 1030
#define DAS_D_KEY 1031
#define DAS_H0_KEY 1032
#define DAS_D0_KEY 1033
#define DAS_HDRSIZE 200000
#define DAS_BUFSIZE 545259520

typedef struct {
  int s0, s1, card;
  int delay, p0, pstep0, p1, pstep1, fstc;
  float p2fstc, fstc_step;
} FftDelayParType;

typedef struct {
  double clock, t_update;
  double pc_time;
  double t_off;
  double delay_step, fstc_scale, nco_tick_time;
  int cycle, seq, cards;
  unsigned char dpcmux, clksel, fftmode, macmode;
  ModelParType par[MAX_SAMPS];
  FftDelayParType fdpar[MAX_SAMPS / 2];
} ModelInfoType;

typedef struct {
  int active, status, scan, scan_off;
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
  unsigned char type, cmd;
  unsigned short seq;
  int flag_num, scan_num;
} EventLogType;

typedef struct {
  int t_off, wt_off, par_off, data_off, data_words;
  int par_words, wordsize;
} RecOffsetType;

typedef struct {
  int off;
  BaseParType base;
  char name[12];
} MacWordType;

typedef struct {
  unsigned char seq, scan;
  int status, recl, seqnum;
  RecOffsetType off;
  MacWordType *mac;
  char *buf;
} RecBufType;

typedef struct {
  int active, status;
  unsigned short events, flags, starts, stops;
  CorrType corr;
  AntennaFlagType flag[MAX_EVENTS][MAX_BANDS];
  EventLogType event[MAX_EVENTS];
  ScanInfoType scaninfo[MAX_SCANS];
  RecOffsetType offset;
} DataInfoType;

typedef struct {
  int flag, rec, seqnum;
  unsigned short flag_seq, newstate;
} DataTabType;

typedef struct {
  int flag, blocksize, maxblocks, cur_block, first_block, cur_rec;
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

constexpr int HDRKEY = 2031;
constexpr int BUFKEY = 2032;

constexpr int BCNTBLKS = 4;
constexpr int MAXBLKS = 12;
constexpr int FFTSAMPS = 800;
constexpr int BLKSAMPS = FFTSAMPS * 32;

typedef struct {
  unsigned int active, status, isempty;
  double pctime, reftime, rectime;
  struct timeval timestamps[MAXBLKS];
  double nanoseconds[MAXBLKS];
  unsigned int flag, curblk, currec, blksize, nbeams;
  int overflow;
} BeamBufferType;

namespace nb = nanobind;
void init_multifrb(nb::module_ m);

#endif
