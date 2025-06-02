#include "multifrb.h"
#include "utilities.h"

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
  /** Attach to header. **/
  m_hdrid = shmget(FRBHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
  if (m_hdrid < 0)
    throw std::runtime_error(
        "Could not obtain header shared memory ID. Exiting...");
  m_hdrptr = (BeamHeaderType *)shmat(m_hdrid, NULL, SHM_RDONLY);
  if ((void *)m_hdrptr == (void *)-1)
    throw std::runtime_error(
        "Could not attach to header shared memory. Exiting...");

  /** Read in all header parameters... **/
  ScanInfoType *scan = &(m_hdrptr->ScanTab[0]);

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
  m_dt = m_hdrptr->corr.daspar.gsb_final_bw *
         m_hdrptr->BeamGenHdr.SampInterval / (m_hdrptr->corr.corrpar.clock);

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
  m_nbeamspernode =
      m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode;

  /** Get beam RA and DEC values. **/
  for (int i = 0; i < m_nbeamspernode; i++) {
    int b = m_beamid * m_nbeamspernode + i;
    m_beamras.push_back(m_hdrptr->BeamGenHdr.BeamSteeringParams.RA[b]);
    m_beamdecs.push_back(m_hdrptr->BeamGenHdr.BeamSteeringParams.DEC[b]);
  }

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

  m_obsdate = safe_strftime("%d/%m/%Y", lt);
  m_obstime = safe_fmt("%02d:%02d:%012.9lf", lt->tm_hour, lt->tm_min, seconds);
}

void MFS::unlink() {
  if (shmdt(m_hdrptr) == -1)
    throw std::runtime_error(
        "Could not detach from header shared memory. Exiting...");
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
  return getslice(beam, t0, t0 + dm2delay(fl(), fh(), dm) + width);
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

      /** Summarise all properties in a dictionary. **/
      .def_prop_ro("header",
                   [](MFS &x) {
                     nb::dict header;
                     header["nf"] = x.nf();
                     header["fh"] = x.fh();
                     header["fl"] = x.fl();
                     header["df"] = x.df();
                     header["bw"] = x.bw();
                     header["dt"] = x.dt();
                     header["ra"] = x.ra();
                     header["dec"] = x.dec();
                     header["nbits"] = x.nbits();
                     header["beamid"] = x.beamid();
                     header["hostid"] = x.hostid();
                     header["nbeams"] = x.nbeams();
                     header["source"] = x.source();
                     header["nstokes"] = x.nstokes();
                     header["flipped"] = x.flipped();
                     header["obsdate"] = x.obsdate();
                     header["obstime"] = x.obstime();
                     header["beamras"] = x.beamras();
                     header["beamdecs"] = x.beamdecs();
                     header["hostname"] = x.hostname();
                     header["beammode"] = x.beammode();
                     header["observer"] = x.observer();
                     header["antspol1"] = x.antspol1();
                     header["antspol2"] = x.antspol2();
                     header["gtaccode"] = x.gtaccode();
                     header["gtactitle"] = x.gtactitle();
                     header["antmaskpol1"] = x.antmaskpol1();
                     header["antmaskpol2"] = x.antmaskpol2();
                     header["npcbaselines"] = x.npcbaselines();
                     header["nbeamspernode"] = x.nbeamspernode();
                     return header;
                   })

      /** PART IV: Shared memory properties. **/
      .def_prop_ro("size", [](MFS &x) { return x.size(); })
      .def_prop_ro("maxblks", [](MFS &x) { return x.maxblks(); })
      .def_prop_ro("blksize", [](MFS &x) { return x.blksize(); })
      .def_prop_ro("blksamps", [](MFS &x) { return x.blksamps(); })
      .def_prop_ro("blktime", [](MFS &x) { return x.blktime(); })
      .def_prop_ro("curtime", [](MFS &x) { return x.curtime(); })
      .def_prop_ro("begtime", [](MFS &x) { return x.begtime(); })
      .def_prop_ro("endtime", [](MFS &x) { return x.endtime(); })
      .def_prop_ro("currec", [](MFS &x) { return x.currec(); })
      .def_prop_ro("curblk", [](MFS &x) { return x.curblk(); })

      /** Class methods. **/
      .def("link", &MFS::link)
      .def("unlink", &MFS::unlink)
      .def("getblk", &MFS::getblk)
      .def("getslice", &MFS::getslice)
      .def("getburst", &MFS::getburst)
      .def("timeofblk", &MFS::timeofblk);
}
