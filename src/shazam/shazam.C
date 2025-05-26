#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "multifrb.h"

namespace nb = nanobind;
using MFS = MultiFRBSHM;

NB_MODULE(internals, m) {
  nb::class_<MFS>(m, "MultiFRBSHM")
      .def(nb::init<>())
      .def_prop_ro("nf", [](MFS &x) { return x.nf(); })
      .def_prop_ro("fh", [](MFS &x) { return x.fh(); })
      .def_prop_ro("fl", [](MFS &x) { return x.fl(); })
      .def_prop_ro("df", [](MFS &x) { return x.df(); })
      .def_prop_ro("bw", [](MFS &x) { return x.bw(); })
      .def_prop_ro("dt", [](MFS &x) { return x.dt(); })
      .def_prop_ro("ra", [](MFS &x) { return x.ra(); })
      .def_prop_ro("dec", [](MFS &x) { return x.dec(); })
      .def_prop_ro("nbits", [](MFS &x) { return x.nbits(); })
      .def_prop_ro("beamid", [](MFS &x) { return x.beamid(); })
      .def_prop_ro("hostid", [](MFS &x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](MFS &x) { return x.nbeams(); })
      .def_prop_ro("nstokes", [](MFS &x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](MFS &x) { return x.flipped(); })
      .def_prop_ro("source", [](MFS &x) { return x.source(); })
      .def_prop_ro("obsdate", [](MFS &x) { return x.obsdate(); })
      .def_prop_ro("obstime", [](MFS &x) { return x.obstime(); })
      .def_prop_ro("hostname", [](MFS &x) { return x.hostname(); })
      .def_prop_ro("beammode", [](MFS &x) { return x.beammode(); })
      .def_prop_ro("observer", [](MFS &x) { return x.observer(); })
      .def_prop_ro("npcbaselines", [](MFS &x) { return x.npcbaselines(); })
      .def_prop_ro("gtaccode", [](MFS &x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](MFS &x) { return x.gtactitle(); })
      .def_prop_ro("nbeamspernode", [](MFS &x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](MFS &x) { return x.beamras(); })
      .def_prop_ro("antmaskpol1", [](MFS &x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](MFS &x) { return x.antmaskpol2(); })
      .def_prop_ro("beamdecs", [](MFS &x) { return x.beamdecs(); })
      .def_prop_ro("antspol1", [](MFS &x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](MFS &x) { return x.antspol2(); })
      .def_prop_ro("header", [](MFS &x) {
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
      });
}
