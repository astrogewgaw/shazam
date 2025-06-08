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

#endif
