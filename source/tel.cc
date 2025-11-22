#include "../include/shazam/tel.h"

#include <stdexcept>

unsigned char* TELRing::ptrtobeam(int beam) {
  if (m_linked) return m_dataptr + blksize() * beam;
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}

unsigned char* TELRing::ptrtoblk(int beam, int blk) {
  if (m_linked) return ptrtobeam(beam) + (blksize() * m_nbeamspernode * (blk % maxblks()));
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}

unsigned char* TELRing::ptrtotime(int beam, double t) {
  if (m_linked) {
    if (t > curtime()) throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
    int blk = (int)std::floor(t / blktime());
    int leftsamps = (int)std::round((t - blk * blktime()) / m_dt);
    return ptrtoblk(beam, blk) + (long)leftsamps * (long)m_nf;
  }
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}

void TELRing::link() {
  if (not m_linked) {
    /** Link the header. **/
    m_header.link();

    /** Transfer some private variables from header instance. **/
    m_hdrid = m_header.m_hdrid;
    m_hdrptr = m_header.m_hdrptr;

    /** Transfer all metadata from the header instance. **/
    m_nf = m_header.m_nf;
    m_fh = m_header.m_fh;
    m_fl = m_header.m_fl;
    m_df = m_header.m_df;
    m_bw = m_header.m_bw;
    m_dt = m_header.m_dt;
    m_ra = m_header.m_ra;
    m_dec = m_header.m_dec;
    m_nbits = m_header.m_nbits;
    m_beamid = m_header.m_beamid;
    m_hostid = m_header.m_hostid;
    m_nbeams = m_header.m_nbeams;
    m_source = m_header.m_source;
    m_nstokes = m_header.m_nstokes;
    m_flipped = m_header.m_flipped;
    m_beamras = m_header.m_beamras;
    m_beamdecs = m_header.m_beamdecs;
    m_hostname = m_header.m_hostname;
    m_beammode = m_header.m_beammode;
    m_observer = m_header.m_observer;
    m_antspol1 = m_header.m_antspol1;
    m_antspol2 = m_header.m_antspol2;
    m_gtaccode = m_header.m_gtaccode;
    m_gtactitle = m_header.m_gtactitle;
    m_antmaskpol1 = m_header.m_antmaskpol1;
    m_antmaskpol2 = m_header.m_antmaskpol2;
    m_npcbaselines = m_header.m_npcbaselines;
    m_nbeamspernode = m_header.m_nbeamspernode;

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

    m_bufid = shmget(ShmKey, curshmsize, SHM_RDONLY);
    if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET TEL SHM ID. ABORT.");
    m_bufptr = (GlobalInfoType*)shmat(m_bufid, NULL, SHM_RDONLY);
    if ((void*)m_bufptr == (void*)-1) throw std::runtime_error("FAILED TO LINK TO TEL SHM. ABORT.");
    m_dataptr = (unsigned char*)m_bufptr;

    /** If everything goes well, update status. **/
    m_linked = true;
  }
}

void TELRing::unlink() {
  if (m_linked) {
    m_header.unlink();
    if (shmdt(m_bufptr) == -1) throw std::runtime_error("FAILED TO UNLINK FROM TEL SHM. ABORT.");
    m_linked = false;
  }
}

std::tuple<unsigned char*, size_t> TELRing::getblk_unsafe(int beam, int blk) {
  unsigned char* ptr = ptrtoblk(beam, blk);
  size_t size = blksamps() * m_nf;
  unsigned char* buffer = new unsigned char[size];
  for (int i = 0; i < blksize(); ++i) buffer[i] = ptr[i];
  return std::make_tuple(buffer, size);
}

std::tuple<unsigned char*, size_t> TELRing::getblk(int beam, int blk) {
  if (m_linked) {
    if (timeofblk(blk) > curtime()) throw std::runtime_error("BLOCK NOT YET WRITTEN. ABORT.");
    return getblk_unsafe(beam, blk);
  }
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}

std::tuple<unsigned char*, size_t> TELRing::getblks_unsafe(int beam, int blk0, int blkN) {
  int nblks = blkN - blk0 + 1;
  size_t size = (size_t)nblks * nf();
  unsigned char* buffer = new unsigned char[nblks * blksamps() * m_nf];
  for (int iblk = 0; iblk < nblks; ++iblk) {
    unsigned char* ptr = ptrtoblk(beam, blk0 + iblk);
    for (int i = iblk * blksize(); i < (iblk + 1) * blksize(); ++i) buffer[i] = ptr[i];
  }
  return std::make_tuple(buffer, size);
}

std::tuple<unsigned char*, size_t> TELRing::getblks(int beam, int blk0, int blkN) {
  if (m_linked) {
    if (timeofblk(blk0) > curtime()) throw std::runtime_error("1ST BLOCK NOT YET WRITTEN. ABORT.");
    if (timeofblk(blkN) > curtime()) throw std::runtime_error("NTH BLOCK NOT YET WRITTEN. ABORT.");
    return getblks_unsafe(beam, blk0, blkN);
  }
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}

std::tuple<unsigned char*, size_t> TELRing::getslice_unsafe(int beam, double tbeg, double tend) {
  size_t begN = (size_t)std::round(tbeg / m_dt);
  size_t endN = (size_t)std::round(tend / m_dt);
  size_t N = endN - begN;
  size_t size = N * m_nf;

  unsigned char* buffer = new unsigned char[size];

  int blk = (int)std::floor(tbeg / blktime());
  unsigned char* ptr = ptrtotime(beam, tbeg);
  unsigned char* endptr = ptrtotime(beam, tend);
  unsigned char* blkptr = ptrtoblk(beam, blk) + blksize();

  for (size_t i = 0;; ++i, ++ptr) {
    if (ptr == blkptr) {
      blk += 1;
      ptr = ptrtoblk(beam, blk);
      blkptr = ptrtoblk(beam, blk) + blksize();
    }
    if (ptr == endptr) break;
    buffer[i] = *ptr;
  }

  return std::make_tuple(buffer, size);
}

std::tuple<unsigned char*, size_t> TELRing::getslice(int beam, double tbeg, double tend) {
  if (m_linked) {
    if ((tbeg > curtime()) || (tend > curtime()))
      throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
    if (curtime() >= ((unsigned int)std::floor(tbeg / blktime()) + maxblks()) * blktime())
      throw std::runtime_error("DATA OVERWRITTEN. ABORT.");
    return getslice_unsafe(beam, tbeg, tend);
  }
  throw std::runtime_error("NO LINK TO TEL SHM. ABORT.");
}
