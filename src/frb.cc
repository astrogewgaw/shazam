#include "../include/shazam/frb.h"

#include <cstring>
#include <stdexcept>
#include <tuple>

namespace shazam {
  constexpr double KDM = 1 / 2.41e-4;

  unsigned char* FRBRing::ptrtobeam(int beam) {
    if (m_opened) return m_dataptr + blksize() * beam;
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  unsigned char* FRBRing::ptrtoblk(int beam, int blk) {
    if (m_opened) return ptrtobeam(beam) + (blksize() * nbeamspernode() * (blk % maxblks()));
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  unsigned char* FRBRing::ptrtotime(int beam, double t) {
    if (m_opened) {
      if (t > curtime()) throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
      int blk = (int)std::floor(t / blktime());
      int leftsamps = (int)std::round((t - blk * blktime()) / dt());
      return ptrtoblk(beam, blk) + (long)leftsamps * (long)nf();
    }
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  void FRBRing::putblk(unsigned char* data, int beam, int blk) {
    if (m_opened) {
      size_t size = blksamps() * nf();
      unsigned char* ptr = ptrtoblk(beam, blk);
      for (int i = 0; i < blksize(); ++i) ptr[i] = data[i];
    }
  }

  void FRBRing::putblks(unsigned char* data, int beam, int blk0, int blkN) {
    if (m_opened) {
      int nblks = blkN - blk0 + 1;
      size_t size = (size_t)nblks * nf();
      for (int iblk = 0; iblk < nblks; ++iblk) {
        unsigned char* ptr = ptrtoblk(beam, blk0 + iblk);
        for (int i = iblk * blksize(); i < (iblk + 1) * blksize(); ++i) ptr[i] = data[i];
      }
    }
  }

  std::tuple<unsigned char*, size_t> FRBRing::getblk(int beam, int blk) {
    if (m_opened) {
      if (timeofblk(blk) > curtime()) throw std::runtime_error("BLOCK NOT YET WRITTEN. ABORT.");

      unsigned char* ptr = ptrtoblk(beam, blk);
      size_t size = blksamps() * nf();
      unsigned char* buffer = new unsigned char[size];
      for (int i = 0; i < blksize(); ++i) buffer[i] = ptr[i];
      return std::make_tuple(buffer, size);
    }
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  std::tuple<unsigned char*, size_t> FRBRing::getblks(int beam, int blk0, int blkN) {
    if (m_opened) {
      if (timeofblk(blk0) > curtime())
        throw std::runtime_error("1ST BLOCK NOT YET WRITTEN. ABORT.");
      if (timeofblk(blkN) > curtime())
        throw std::runtime_error("NTH BLOCK NOT YET WRITTEN. ABORT.");

      int nblks = blkN - blk0 + 1;
      size_t size = (size_t)nblks * nf();
      unsigned char* buffer = new unsigned char[nblks * blksamps() * nf()];
      for (int iblk = 0; iblk < nblks; ++iblk) {
        unsigned char* ptr = ptrtoblk(beam, blk0 + iblk);
        for (int i = iblk * blksize(); i < (iblk + 1) * blksize(); ++i) buffer[i] = ptr[i];
      }
      return std::make_tuple(buffer, size);
    }
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  std::tuple<unsigned char*, size_t> FRBRing::getslice(int beam, double tbeg, double tend) {
    if (m_opened) {
      if ((tbeg > curtime()) || (tend > curtime()))
        throw std::runtime_error("DATA NOT YET WRITTEN. ABORT.");
      if (curtime() >= ((unsigned int)std::floor(tbeg / blktime()) + maxblks()) * blktime())
        throw std::runtime_error("DATA OVERWRITTEN. ABORT.");

      size_t begN = (size_t)std::round(tbeg / dt());
      size_t endN = (size_t)std::round(tend / dt());
      size_t N = endN - begN;
      size_t size = N * nf();

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
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }

  std::tuple<unsigned char*, size_t> FRBRing::getburst(int beam, double t0, double dm,
                                                       double width) {
    if (m_opened) {
      double delay = KDM * dm * (std::pow(fl(), -2) - std::pow(fh(), -2));
      double tend = t0 + width + delay;
      double tbeg = t0 - width;
      return getslice(beam, tbeg, tend);
    }
    throw std::runtime_error("FRB SHM NOT OPEN. ABORT.");
  }
}  // namespace shazam
