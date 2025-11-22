#include "../include/shazam/hdr.h"

#include <stdexcept>

void Header::link() {
  if (not m_linked) {
    /** Attach to header. **/
    m_hdrid = shmget(MULTIHDRKEY, sizeof(BeamHeaderType), SHM_RDONLY);
    if (m_hdrid < 0) throw std::runtime_error("UNABLE TO GET HDR SHM ID. ABORT.");
    m_hdrptr = (BeamHeaderType*)shmat(m_hdrid, NULL, SHM_RDONLY);
    if ((void*)m_hdrptr == (void*)-1) throw std::runtime_error("FAILED TO LINK TO HDR SHM. ABORT");

    /** Read in all header parameters... **/
    ScanInfoType* scan = &(m_hdrptr->ScanTab[0]);

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
    m_dt = m_hdrptr->corr.daspar.gsb_final_bw * m_hdrptr->BeamGenHdr.SampInterval
           / (m_hdrptr->corr.corrpar.clock);

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
    m_nbeamspernode = m_hdrptr->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode;

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

void Header::unlink() {
  if (m_linked) {
    if (shmdt(m_hdrptr) == -1) throw std::runtime_error("FAILED TO UNLINK FROM HDR SHM. ABORT.");
    m_linked = false;
  }
}
