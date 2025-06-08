#include "multifrb.h"

using MFS = MultiFRBSHM;

unsigned char *MFS::ptrtobeam(int beam) { return m_dataptr + blksize() * beam; }

unsigned char *MFS::ptrtoblk(int beam, int blk) {
  return ptrtobeam(beam) + (blksize() * nbeamspernode() * (blk % maxblks()));
}

unsigned char *MFS::ptrtotime(int beam, double t) {
  if (t > curtime())
    throw std::runtime_error("Data not yet written. Exiting...");
  int blk = (int)std::floor(t / blktime());
  int leftsamps = (int)std::round((t - blk * blktime()) / dt());
  return ptrtoblk(beam, blk) + (long)leftsamps * (long)nf();
}

void MFS::link() {
  m_hdrid = m_header.m_hdrid;
  m_hdrptr = m_header.m_hdrptr;

  m_nf = m_header.nf();
  m_fh = m_header.fh();
  m_fl = m_header.fl();
  m_df = m_header.df();
  m_bw = m_header.bw();
  m_dt = m_header.dt();
  m_ra = m_header.ra();
  m_dec = m_header.dec();
  m_nbits = m_header.nbits();
  m_beamid = m_header.beamid();
  m_hostid = m_header.hostid();
  m_nbeams = m_header.nbeams();
  m_source = m_header.source();
  m_nstokes = m_header.nstokes();
  m_flipped = m_header.flipped();
  m_beamras = m_header.beamras();
  m_beamdecs = m_header.beamdecs();
  m_hostname = m_header.hostname();
  m_beammode = m_header.beammode();
  m_observer = m_header.observer();
  m_antspol1 = m_header.antspol1();
  m_antspol2 = m_header.antspol2();
  m_gtaccode = m_header.gtaccode();
  m_gtactitle = m_header.gtactitle();
  m_antmaskpol1 = m_header.antmaskpol1();
  m_antmaskpol2 = m_header.antmaskpol2();
  m_npcbaselines = m_header.npcbaselines();
  m_nbeamspernode = m_header.nbeamspernode();

  /** Calculate size of buffer. **/
  long BLKSIZE = (long)FRBBLKSAMPS * (long)m_nf;
  long BUFSIZE = BLKSIZE * (long)FRBMAXBLKS * (long)m_nbeamspernode;
  long FRBSHMSIZE = sizeof(BeamBufferType) + BUFSIZE;

  /** Attach to buffer. **/
  m_bufid = shmget(FRBBUFKEY, FRBSHMSIZE, SHM_RDONLY);
  if (m_bufid < 0)
    throw std::runtime_error(
        "Could not obtain buffer shared memory ID. Exiting...");
  m_bufptr = (BeamBufferType *)shmat(m_bufid, NULL, SHM_RDONLY);
  if ((void *)m_bufptr == (void *)-1)
    throw std::runtime_error(
        "Could not attach to buffer shared memory. Exiting...");
  m_dataptr = ((unsigned char *)m_bufptr) + sizeof(BeamBufferType);

  /** Get index of current record and block being written. **/
  m_curblk = m_bufptr->curblk;
  m_currec = (m_bufptr->empty) ? m_bufptr->currec
                               : (m_bufptr->currec - 1) % FRBMAXBLKS;

  /** Get date and time of observation. **/
  struct tm *lt;
  double seconds;
  double nanoseconds = 0;
  struct timeval timestamp;

  memcpy(&timestamp, &(m_bufptr->timestamps[m_currec]), sizeof(struct timeval));
  memcpy(&nanoseconds, &(m_bufptr->nanoseconds[m_currec]), sizeof(nanoseconds));

  lt = localtime(&timestamp.tv_sec);
  seconds = lt->tm_sec + (timestamp.tv_usec) / 1e6 + (nanoseconds) / 1e9;

  m_obsdate = sstrftime("%d/%m/%Y", lt);
  m_obstime = sfmt("%02d:%02d:%012.9lf", lt->tm_hour, lt->tm_min, seconds);
}

void MFS::unlink() {
  m_header.unlink();
  if (shmdt(m_bufptr) == -1)
    throw std::runtime_error(
        "Could not detach from buffer shared memory. Exiting...");
}

Array MFS::getblk(int beam, int blk) {
  if (timeofblk(blk) > curtime())
    throw std::runtime_error("Block not yet written. Exiting...");
  unsigned char *ptr = ptrtoblk(beam, blk);
  unsigned char *buffer = new unsigned char[blksamps() * m_nf];
  for (int i = 0; i < blksize(); ++i) buffer[i] = ptr[i];
  return Array(buffer, {(size_t)blksamps(), (size_t)nf()},
               nb::capsule(buffer, [](void *p) noexcept {
                 delete[] (unsigned char *)p;
               }));
}

Array MFS::getblks(int beam, int blk0, int blkN) {
  if (timeofblk(blk0) > curtime())
    throw std::runtime_error("First block not yet written. Exiting...");
  if (timeofblk(blkN) > curtime())
    throw std::runtime_error("Last block not yet written. Exiting...");

  int nblks = blkN - blk0 + 1;
  unsigned char *buffer = new unsigned char[nblks * blksamps() * m_nf];
  for (int iblk = 0; iblk < nblks; ++iblk) {
    unsigned char *ptr = ptrtoblk(beam, blk0 + iblk);
    for (int i = iblk * blksize(); i < (iblk + 1) * blksize(); ++i)
      buffer[i] = ptr[i];
  }

  return Array(buffer, {(size_t)(nblks * blksamps()), (size_t)nf()},
               nb::capsule(buffer, [](void *p) noexcept {
                 delete[] (unsigned char *)p;
               }));
}

Array MFS::getslice(int beam, double tbeg, double tend) {
  if ((tbeg > curtime()) || (tend > curtime()))
    throw std::runtime_error("Data has not yet been written. Exiting...");
  if (curtime() >=
      ((unsigned int)std::floor(tbeg / blktime()) + 12) * blktime())
    throw std::runtime_error("Data has been overwritten. Exiting...");

  size_t begN = (size_t)std::round(tbeg / m_dt);
  size_t endN = (size_t)std::round(tend / m_dt);
  size_t N = endN - begN;

  unsigned char *buffer = new unsigned char[N * m_nf];

  unsigned char *ptr = ptrtotime(beam, tbeg);
  unsigned char *endptr = ptrtotime(beam, tend);
  int blk = (int)std::floor(tbeg / blktime());
  unsigned char *blkptr = ptrtoblk(beam, blk) + blksize();

  for (size_t i = 0;; ++i, ++ptr) {
    if (ptr == blkptr) {
      blk += 1;
      ptr = ptrtoblk(beam, blk);
      blkptr = ptrtoblk(beam, blk) + blksize();
    }
    if (ptr == endptr) break;
    buffer[i] = *ptr;
  }

  return Array(buffer, {N, (size_t)m_nf},
               nb::capsule(buffer, [](void *p) noexcept {
                 delete[] (unsigned char *)p;
               }));
}

Array MFS::getburst(int beam, double t0, double dm, double width) {
  double tbeg = t0 - width;
  double tend = t0 + width + dm2delay(fl(), fh(), dm);
  return getslice(beam, tbeg, tend);
}

void initmultifrb(nb::module_ m) {
  nb::class_<MFS>(m, "MultiFRBSHM")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](MFS &x) { return x.nf(); })
      .def_prop_ro("fh", [](MFS &x) { return x.fh(); })
      .def_prop_ro("fl", [](MFS &x) { return x.fl(); })
      .def_prop_ro("df", [](MFS &x) { return x.df(); })
      .def_prop_ro("bw", [](MFS &x) { return x.bw(); })
      .def_prop_ro("dt", [](MFS &x) { return x.dt(); })
      .def_prop_ro("nbits", [](MFS &x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](MFS &x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](MFS &x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](MFS &x) { return x.ra(); })
      .def_prop_ro("dec", [](MFS &x) { return x.dec(); })
      .def_prop_ro("source", [](MFS &x) { return x.source(); })
      .def_prop_ro("obsdate", [](MFS &x) { return x.obsdate(); })
      .def_prop_ro("obstime", [](MFS &x) { return x.obstime(); })
      .def_prop_ro("beammode", [](MFS &x) { return x.beammode(); })
      .def_prop_ro("observer", [](MFS &x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](MFS &x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](MFS &x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](MFS &x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](MFS &x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](MFS &x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](MFS &x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](MFS &x) { return x.beamid(); })
      .def_prop_ro("hostid", [](MFS &x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](MFS &x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](MFS &x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](MFS &x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](MFS &x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](MFS &x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](MFS &x) { return x.beamdecs(); })

      .def_prop_ro("header", [](MFS &x) { return x.header(); })
      .def_prop_ro("attrs", [](MFS &x) { return x.header().asdict(); })

      /** PART IV: Shared memory properties. **/
      .def_prop_ro("size", [](MFS &x) { return x.size(); })
      .def_prop_ro("empty", [](MFS &x) { return x.empty(); })
      .def_prop_ro("status", [](MFS &x) { return x.status(); })
      .def_prop_ro("active", [](MFS &x) { return x.active(); })
      .def_prop_ro("maxblks", [](MFS &x) { return x.maxblks(); })
      .def_prop_ro("blksize", [](MFS &x) { return x.blksize(); })
      .def_prop_ro("blksamps", [](MFS &x) { return x.blksamps(); })
      .def_prop_ro("blktime", [](MFS &x) { return x.blktime(); })
      .def_prop_ro("curtime", [](MFS &x) { return x.curtime(); })
      .def_prop_ro("begtime", [](MFS &x) { return x.begtime(); })
      .def_prop_ro("endtime", [](MFS &x) { return x.endtime(); })
      .def_prop_ro("currec", [](MFS &x) { return x.currec(); })
      .def_prop_ro("curblk", [](MFS &x) { return x.curblk(); })

      /** Public methods. **/
      .def("link", &MFS::link)
      .def("unlink", &MFS::unlink)
      .def("timeofblk", &MFS::timeofblk, "blk"_a)
      .def("getblk", &MFS::getblk, "beam"_a, "blk"_a)
      .def("getblks", &MFS::getblks, "beam"_a, "blk0"_a, "blkN"_a)
      .def("getslice", &MFS::getslice, "beam"_a, "tbeg"_a, "tend"_a)
      .def("getburst", &MFS::getburst, "beam"_a, "t0"_a, "dm"_a, "width"_a);
}
