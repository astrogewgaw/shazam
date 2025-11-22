#include "../include/shazam/frb.h"

#include <cstring>
#include <stdexcept>
#include <tuple>

constexpr double KDM = 1 / 2.41e-4;

unsigned char* FRBRing::ptrtobeam(int beam) {
  if (m_linked) return m_dataptr + blksize() * beam;
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

unsigned char* FRBRing::ptrtoblk(int beam, int blk) {
  if (m_linked) return ptrtobeam(beam) + (blksize() * m_nbeamspernode * (blk % maxblks()));
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

unsigned char* FRBRing::ptrtotime(int beam, double t) {
  if (m_linked) {
    if (t > curtime()) throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
    int blk = (int)std::floor(t / blktime());
    int leftsamps = (int)std::round((t - blk * blktime()) / m_dt);
    return ptrtoblk(beam, blk) + (long)leftsamps * (long)m_nf;
  }
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

void FRBRing::link() {
  if (not m_linked) {
    /** Link the header. **/
    m_hdr.link();

    /** Transfer some private variables from header instance. **/
    m_hdrid = m_hdr.m_hdrid;
    m_hdrptr = m_hdr.m_hdrptr;

    /** Transfer all metadata from the header instance. **/
    m_nf = m_hdr.m_nf;
    m_fh = m_hdr.m_fh;
    m_fl = m_hdr.m_fl;
    m_df = m_hdr.m_df;
    m_bw = m_hdr.m_bw;
    m_dt = m_hdr.m_dt;
    m_ra = m_hdr.m_ra;
    m_dec = m_hdr.m_dec;
    m_nbits = m_hdr.m_nbits;
    m_beamid = m_hdr.m_beamid;
    m_hostid = m_hdr.m_hostid;
    m_nbeams = m_hdr.m_nbeams;
    m_source = m_hdr.m_source;
    m_nstokes = m_hdr.m_nstokes;
    m_flipped = m_hdr.m_flipped;
    m_beamras = m_hdr.m_beamras;
    m_beamdecs = m_hdr.m_beamdecs;
    m_hostname = m_hdr.m_hostname;
    m_beammode = m_hdr.m_beammode;
    m_observer = m_hdr.m_observer;
    m_antspol1 = m_hdr.m_antspol1;
    m_antspol2 = m_hdr.m_antspol2;
    m_gtaccode = m_hdr.m_gtaccode;
    m_gtactitle = m_hdr.m_gtactitle;
    m_antmaskpol1 = m_hdr.m_antmaskpol1;
    m_antmaskpol2 = m_hdr.m_antmaskpol2;
    m_npcbaselines = m_hdr.m_npcbaselines;
    m_nbeamspernode = m_hdr.m_nbeamspernode;

    /** Calculate size of buffer. **/
    long BLKSIZE = (long)FRBBLKSAMPS * (long)m_nf;
    long BUFSIZE = BLKSIZE * (long)FRBMAXBLKS * (long)m_nbeamspernode;
    long FRBSHMSIZE = sizeof(BeamBufferType) + BUFSIZE;

    /** Attach to buffer. **/
    m_bufid = shmget(FRBBUFKEY, FRBSHMSIZE, SHM_RDONLY);
    if (m_bufid < 0) throw std::runtime_error("UNABLE TO GET FRB SHM ID. ABORT.");
    m_bufptr = (BeamBufferType*)shmat(m_bufid, NULL, SHM_RDONLY);
    if ((void*)m_bufptr == (void*)-1) throw std::runtime_error("FAILED TO LINK TO FRB SHM. ABORT.");
    m_dataptr = ((unsigned char*)m_bufptr) + sizeof(BeamBufferType);

    /** If everything goes well, update status. **/
    m_linked = true;
  }
}

void FRBRing::unlink() {
  if (m_linked) {
    m_hdr.unlink();
    if (shmdt(m_bufptr) == -1) throw std::runtime_error("FAILED TO UNLINK FROM FRB SHM. ABORT.");
    m_linked = false;
  }
}

std::tuple<unsigned char*, size_t> FRBRing::getblk_unsafe(int beam, int blk) {
  unsigned char* ptr = ptrtoblk(beam, blk);
  size_t size = blksamps() * m_nf;
  unsigned char* buffer = new unsigned char[size];
  for (int i = 0; i < blksize(); ++i) buffer[i] = ptr[i];
  return std::make_tuple(buffer, size);
}

std::tuple<unsigned char*, size_t> FRBRing::getblk(int beam, int blk) {
  if (m_linked) {
    if (timeofblk(blk) > curtime()) throw std::runtime_error("BLOCK NOT YET WRITTEN. ABORT.");
    return getblk_unsafe(beam, blk);
  }
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

std::tuple<unsigned char*, size_t> FRBRing::getblks_unsafe(int beam, int blk0, int blkN) {
  int nblks = blkN - blk0 + 1;
  size_t size = (size_t)nblks * nf();
  unsigned char* buffer = new unsigned char[nblks * blksamps() * m_nf];
  for (int iblk = 0; iblk < nblks; ++iblk) {
    unsigned char* ptr = ptrtoblk(beam, blk0 + iblk);
    for (int i = iblk * blksize(); i < (iblk + 1) * blksize(); ++i) buffer[i] = ptr[i];
  }
  return std::make_tuple(buffer, size);
}

std::tuple<unsigned char*, size_t> FRBRing::getblks(int beam, int blk0, int blkN) {
  if (m_linked) {
    if (timeofblk(blk0) > curtime()) throw std::runtime_error("1ST BLOCK NOT YET WRITTEN. ABORT.");
    if (timeofblk(blkN) > curtime()) throw std::runtime_error("NTH BLOCK NOT YET WRITTEN. ABORT.");
    return getblks_unsafe(beam, blk0, blkN);
  }
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

std::tuple<unsigned char*, size_t> FRBRing::getslice_unsafe(int beam, double tbeg, double tend) {
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

std::tuple<unsigned char*, size_t> FRBRing::getslice(int beam, double tbeg, double tend) {
  if (m_linked) {
    if ((tbeg > curtime()) || (tend > curtime()))
      throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
    if (curtime() >= ((unsigned int)std::floor(tbeg / blktime()) + maxblks()) * blktime())
      throw std::runtime_error("DATA OVERWRITTEN. ABORT.");
    return getslice_unsafe(beam, tbeg, tend);
  }
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}

std::tuple<unsigned char*, size_t> FRBRing::getburst(int beam, double t0, double dm, double width) {
  if (m_linked) {
    double delay = KDM * dm * (std::pow(m_fl, -2) - std::pow(m_fh, -2));
    double tend = t0 + width + delay;
    double tbeg = t0 - width;
    return getslice(beam, tbeg, tend);
  }
  throw std::runtime_error("NO LINK TO FRB SHM. ABORT.");
}
