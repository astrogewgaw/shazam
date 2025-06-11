#include "multihdr.h"

using MH = MultiHeader;

void MH::link() {
  if (not m_linked) {
    /** Attach to header. **/
    m_hdrid = shmget(MULTIHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
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

    /** If everything goes well, update status. **/
    m_linked = true;
  }
}

void MH::unlink() {
  if (m_linked) {
    if (shmdt(m_hdrptr) == -1)
      throw std::runtime_error(
          "Could not detach from header shared memory. Exiting...");
    m_linked = false;
  }
}

nb::dict MH::asdict() {
  nb::dict header;
  header["nf"] = nf();
  header["fh"] = fh();
  header["fl"] = fl();
  header["df"] = df();
  header["bw"] = bw();
  header["dt"] = dt();
  header["ra"] = ra();
  header["dec"] = dec();
  header["nbits"] = nbits();
  header["beamid"] = beamid();
  header["hostid"] = hostid();
  header["nbeams"] = nbeams();
  header["source"] = source();
  header["nstokes"] = nstokes();
  header["flipped"] = flipped();
  header["beamras"] = beamras();
  header["beamdecs"] = beamdecs();
  header["hostname"] = hostname();
  header["beammode"] = beammode();
  header["observer"] = observer();
  header["antspol1"] = antspol1();
  header["antspol2"] = antspol2();
  header["gtaccode"] = gtaccode();
  header["gtactitle"] = gtactitle();
  header["antmaskpol1"] = antmaskpol1();
  header["antmaskpol2"] = antmaskpol2();
  header["npcbaselines"] = npcbaselines();
  header["nbeamspernode"] = nbeamspernode();
  return header;
}

void initmultihdr(nb::module_ m) {
  nb::class_<MH>(m, "MultiHeader")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](MH &x) { return x.nf(); })
      .def_prop_ro("fh", [](MH &x) { return x.fh(); })
      .def_prop_ro("fl", [](MH &x) { return x.fl(); })
      .def_prop_ro("df", [](MH &x) { return x.df(); })
      .def_prop_ro("bw", [](MH &x) { return x.bw(); })
      .def_prop_ro("dt", [](MH &x) { return x.dt(); })
      .def_prop_ro("nbits", [](MH &x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](MH &x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](MH &x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](MH &x) { return x.ra(); })
      .def_prop_ro("dec", [](MH &x) { return x.dec(); })
      .def_prop_ro("source", [](MH &x) { return x.source(); })
      .def_prop_ro("beammode", [](MH &x) { return x.beammode(); })
      .def_prop_ro("observer", [](MH &x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](MH &x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](MH &x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](MH &x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](MH &x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](MH &x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](MH &x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](MH &x) { return x.beamid(); })
      .def_prop_ro("hostid", [](MH &x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](MH &x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](MH &x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](MH &x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](MH &x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](MH &x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](MH &x) { return x.beamdecs(); })

      /** PART IV: Shared memory properties. **/
      .def_prop_ro("linked", [](MH &x) { return x.linked(); })

      /** Dunder methods. **/
      .def("__exit__", [](MH &x, nb::args _) { x.unlink(); })
      .def("__enter__",
           [](MH &x) {
             x.link();
             return x;
           })

      /** Public methods. **/
      .def("link", &MH::link)
      .def("unlink", &MH::unlink)
      .def("asdict", &MH::asdict);
}
