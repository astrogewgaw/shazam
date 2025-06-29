#include "multitel.h"

using MTS = MultiTELSHM;

unsigned char *MTS::ptrtobeam(int beam) {
  if (m_linked) return m_dataptr + blksize() * beam;
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

unsigned char *MTS::ptrtoblk(int beam, int blk) {
  if (m_linked)
    return ptrtobeam(beam) + (blksize() * nbeamspernode() * (blk % maxblks()));
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

unsigned char *MTS::ptrtotime(int beam, double t) {
  if (m_linked) {
    if (t > curtime())
      throw std::runtime_error("Data not yet written. Exiting...");
    int blk = (int)std::floor(t / blktime());
    int leftsamps = (int)std::round((t - blk * blktime()) / dt());
    return ptrtoblk(beam, blk) + (long)leftsamps * (long)nf();
  }
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

void MTS::link() {
  if (not m_linked) {
    /** Link the header. **/
    m_header.link();

    /** Transfer some private variables from header instance. **/
    m_hdrid = m_header.m_hdrid;
    m_hdrptr = m_header.m_hdrptr;

    /** Transfer all metadata from the header instance. **/
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

    int extrabuf = 64;
    int cursamps = 32 * 25;
    long int curtotalwords = cursamps * nf();
    long int currecsize = curtotalwords * WordSize / 2;

    long curshmdatasize =
        (MaxRecs + 1) * currecsize * nbeamspernode() + extrabuf;
    curshmdatasize = curshmdatasize / PageSize + 1;
    curshmdatasize = curshmdatasize * PageSize;

    int shmdataoff = sizeof(GlobalInfoType) + extrabuf;
    shmdataoff = shmdataoff / PageSize + 1;
    shmdataoff = shmdataoff * PageSize;

    long curshmsize = curshmdatasize + shmdataoff;

    m_bufid = shmget(ShmKey, curshmsize, SHM_RDONLY);
    if (m_bufid < 0)
      throw std::runtime_error(
          "Could not obtain buffer shared memory ID. Exiting...");
    m_bufptr = (GlobalInfoType *)shmat(m_bufid, NULL, SHM_RDONLY);
    if ((void *)m_bufptr == (void *)-1)
      throw std::runtime_error(
          "Could not attach to buffer shared memory. Exiting...");
    m_dataptr = (unsigned char *)m_bufptr;

    /** If everything goes well, update status. **/
    m_linked = true;
  }
}

void MTS::unlink() {
  if (m_linked) {
    m_header.unlink();
    if (shmdt(m_bufptr) == -1)
      throw std::runtime_error(
          "Could not detach from buffer shared memory. Exiting...");
    m_linked = false;
  }
}

Array MTS::getblk(int beam, int blk) {
  if (m_linked) {
    if (timeofblk(blk) > curtime())
      throw std::runtime_error("Block not yet written. Exiting...");
    unsigned char *ptr = ptrtoblk(beam, blk);
    unsigned char *buffer = new unsigned char[blksamps() * nf()];
    for (int i = 0; i < blksize(); ++i) buffer[i] = ptr[i];
    return Array(buffer, {(size_t)blksamps(), (size_t)nf()},
                 nb::capsule(buffer, [](void *p) noexcept {
                   delete[] (unsigned char *)p;
                 }));
  }
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

Array MTS::getblks(int beam, int blk0, int blkN) {
  if (m_linked) {
    if (timeofblk(blk0) > curtime())
      throw std::runtime_error("First block not yet written. Exiting...");
    if (timeofblk(blkN) > curtime())
      throw std::runtime_error("Last block not yet written. Exiting...");

    int nblks = blkN - blk0 + 1;
    unsigned char *buffer = new unsigned char[nblks * blksamps() * nf()];
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
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

Array MTS::getslice(int beam, double tbeg, double tend) {
  if (m_linked) {
    if ((tbeg > curtime()) || (tend > curtime()))
      throw std::runtime_error("Data has not yet been written. Exiting...");
    if (curtime() >=
        ((unsigned int)std::floor(tbeg / blktime()) + maxblks()) * blktime())
      throw std::runtime_error("Data has been overwritten. Exiting...");

    size_t begN = (size_t)std::round(tbeg / dt());
    size_t endN = (size_t)std::round(tend / dt());
    size_t N = endN - begN;

    unsigned char *buffer = new unsigned char[N * nf()];

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

    return Array(buffer, {N, (size_t)nf()},
                 nb::capsule(buffer, [](void *p) noexcept {
                   delete[] (unsigned char *)p;
                 }));
  }
  throw std::runtime_error("Not linked to shared memory. Exiting...");
}

void initmultitel(nb::module_ m) {
  nb::class_<MTS>(m, "MultiTELSHM")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](MTS &x) { return x.nf(); })
      .def_prop_ro("fh", [](MTS &x) { return x.fh(); })
      .def_prop_ro("fl", [](MTS &x) { return x.fl(); })
      .def_prop_ro("df", [](MTS &x) { return x.df(); })
      .def_prop_ro("bw", [](MTS &x) { return x.bw(); })
      .def_prop_ro("dt", [](MTS &x) { return x.dt(); })
      .def_prop_ro("nbits", [](MTS &x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](MTS &x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](MTS &x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](MTS &x) { return x.ra(); })
      .def_prop_ro("dec", [](MTS &x) { return x.dec(); })
      .def_prop_ro("source", [](MTS &x) { return x.source(); })
      .def_prop_ro("obsdate", [](MTS &x) { return x.obsdate(); })
      .def_prop_ro("obstime", [](MTS &x) { return x.obstime(); })
      .def_prop_ro("beammode", [](MTS &x) { return x.beammode(); })
      .def_prop_ro("observer", [](MTS &x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](MTS &x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](MTS &x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](MTS &x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](MTS &x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](MTS &x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](MTS &x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](MTS &x) { return x.beamid(); })
      .def_prop_ro("hostid", [](MTS &x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](MTS &x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](MTS &x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](MTS &x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](MTS &x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](MTS &x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](MTS &x) { return x.beamdecs(); })

      /** Summarise all properties as a dictionary. **/
      .def_prop_ro("header", [](MTS &x) { return x.header().asdict(); })

      /** PART IV: Shared memory properties. **/
      .def_prop_ro("linked", [](MTS &x) { return x.linked(); })
      .def_prop_ro("maxblks", [](MTS &x) { return x.maxblks(); })
      .def_prop_ro("blksize", [](MTS &x) { return x.blksize(); })
      .def_prop_ro("blksamps", [](MTS &x) { return x.blksamps(); })
      .def_prop_ro("blktime", [](MTS &x) { return x.blktime(); })
      .def_prop_ro("curtime", [](MTS &x) { return x.curtime(); })
      .def_prop_ro("begtime", [](MTS &x) { return x.begtime(); })
      .def_prop_ro("endtime", [](MTS &x) { return x.endtime(); })
      .def_prop_ro("currec", [](MTS &x) { return x.currec(); })
      .def_prop_ro("curblk", [](MTS &x) { return x.curblk(); })

      /** Dunder methods. **/
      .def("__exit__", [](MTS &x, nb::args _) { x.unlink(); })
      .def("__enter__",
           [](MTS &x) {
             x.link();
             return x;
           })

      /** Public methods. **/
      .def("link", &MTS::link)
      .def("unlink", &MTS::unlink)
      .def("timeofblk", &MTS::timeofblk, "blk"_a)
      .def("getblk", &MTS::getblk, "beam"_a, "blk"_a)
      .def("getblks", &MTS::getblks, "beam"_a, "blk0"_a, "blkN"_a)
      .def("getslice", &MTS::getslice, "beam"_a, "tbeg"_a, "tend"_a);
}
