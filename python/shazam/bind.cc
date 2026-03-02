#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/chrono.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <stdexcept>

#include "../../include/shazam/frb.h"
#include "../../include/shazam/hdr.h"
#include "../../include/shazam/tel.h"

using namespace shazam;
namespace nb = nanobind;
using namespace nb::literals;
using Array = nb::ndarray<nb::numpy, unsigned char, nb::ndim<2>>;

NB_MODULE(core, m) {
  nb::class_<Header>(m, "Header")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      .def_prop_ro("mode",
                   [](Header& x) {
                     switch (x.mode()) {
                       default:
                       case READ:
                         return "r";
                       case WRITE:
                         return "w";
                     }
                   })
      .def_prop_ro("opened", [](Header& x) { return x.opened(); })

      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](Header& x) { return x.nf(); })
      .def_prop_ro("fh", [](Header& x) { return x.fh(); })
      .def_prop_ro("fl", [](Header& x) { return x.fl(); })
      .def_prop_ro("df", [](Header& x) { return x.df(); })
      .def_prop_ro("bw", [](Header& x) { return x.bw(); })
      .def_prop_ro("dt", [](Header& x) { return x.dt(); })
      .def_prop_ro("nbits", [](Header& x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](Header& x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](Header& x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](Header& x) { return x.ra(); })
      .def_prop_ro("dec", [](Header& x) { return x.dec(); })
      .def_prop_ro("source", [](Header& x) { return x.source(); })
      .def_prop_ro("beammode", [](Header& x) { return x.beammode(); })
      .def_prop_ro("observer", [](Header& x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](Header& x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](Header& x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](Header& x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](Header& x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](Header& x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](Header& x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](Header& x) { return x.beamid(); })
      .def_prop_ro("hostid", [](Header& x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](Header& x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](Header& x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](Header& x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](Header& x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](Header& x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](Header& x) { return x.beamdecs(); })

      /** Public methods. **/
      .def("open",
           [](Header& x, std::string mode) {
             if (mode == "r") {
               x.open(READ);
             } else if ((mode == "w") || (mode == "rw")) {
               x.open(WRITE);
             } else {
               throw std::runtime_error("MODE DOESN'T EXIST. ABORT.");
             }
             return x;
           })
      .def("copy", &Header::copy)
      .def("close", &Header::close)
      .def("update", &Header::update)
      .def("asdict", [](Header& x) {
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

  nb::class_<TELRing>(m, "TELRing")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      .def_prop_ro("mode",
                   [](TELRing& x) {
                     switch (x.mode()) {
                       default:
                       case READ:
                         return "r";
                       case WRITE:
                         return "w";
                     }
                   })
      .def_prop_ro("opened", [](TELRing& x) { return x.opened(); })

      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](TELRing& x) { return x.nf(); })
      .def_prop_ro("fh", [](TELRing& x) { return x.fh(); })
      .def_prop_ro("fl", [](TELRing& x) { return x.fl(); })
      .def_prop_ro("df", [](TELRing& x) { return x.df(); })
      .def_prop_ro("bw", [](TELRing& x) { return x.bw(); })
      .def_prop_ro("dt", [](TELRing& x) { return x.dt(); })
      .def_prop_ro("nbits", [](TELRing& x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](TELRing& x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](TELRing& x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](TELRing& x) { return x.ra(); })
      .def_prop_ro("dec", [](TELRing& x) { return x.dec(); })
      .def_prop_ro("source", [](TELRing& x) { return x.source(); })
      .def_prop_ro("beammode", [](TELRing& x) { return x.beammode(); })
      .def_prop_ro("observer", [](TELRing& x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](TELRing& x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](TELRing& x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](TELRing& x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](TELRing& x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](TELRing& x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](TELRing& x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](TELRing& x) { return x.beamid(); })
      .def_prop_ro("hostid", [](TELRing& x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](TELRing& x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](TELRing& x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](TELRing& x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](TELRing& x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](TELRing& x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](TELRing& x) { return x.beamdecs(); })

      /** Summarise all properties as a dictionary. **/
      .def("header",
           [](TELRing& x) {
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
      .def_prop_ro("maxblks", [](TELRing& x) { return x.maxblks(); })
      .def_prop_ro("blksamps", [](TELRing& x) { return x.blksamps(); })

      .def_prop_ro("acqover", [](TELRing& x) { return x.acqover(); })
      .def_prop_ro("gpsok", [](TELRing& x) { return x.gpsok(); })
      .def_prop_ro("acqok", [](TELRing& x) { return x.acqok(); })

      .def_prop_ro("marked", [](TELRing& x) { return x.marked(); })
      .def_prop_ro("dataok", [](TELRing& x) { return x.dataok(); })
      .def_prop_ro("blkok", [](TELRing& x) { return x.blkok(); })
      .def_prop_ro("timeok", [](TELRing& x) { return x.timeok(); })
      .def_prop_ro("noinit", [](TELRing& x) { return x.noinit(); })

      .def_prop_ro("curblk", [](TELRing& x) { return x.curblk(); })
      .def_prop_ro("currec", [](TELRing& x) { return x.currec(); })
      .def_prop_ro("begblk", [](TELRing& x) { return x.begblk(); })
      .def_prop_ro("endblk", [](TELRing& x) { return x.endblk(); })

      .def_prop_ro("blksize", [](TELRing& x) { return x.blksize(); })
      .def_prop_ro("size", [](TELRing& x) { return x.size(); })

      .def_prop_ro("blktime", [](TELRing& x) { return x.blktime(); })
      .def_prop_ro("curtime", [](TELRing& x) { return x.curtime(); })
      .def_prop_ro("begtime", [](TELRing& x) { return x.begtime(); })
      .def_prop_ro("endtime", [](TELRing& x) { return x.endtime(); })
      .def("timeofblk", &TELRing::timeofblk, "blk"_a)
      .def_prop_ro("timestamps", [](FRBRing& x) { return x.timestamps(); })

      /** Public methods. **/
      .def("open",
           [](TELRing& x, std::string mode) {
             if (mode == "r") {
               x.open(READ);
             } else if ((mode == "w") || (mode == "rw")) {
               x.open(WRITE);
             } else {
               throw std::runtime_error("MODE DOESN'T EXIST. ABORT.");
             }
             return x;
           })
      .def("close", &TELRing::close)
      .def("update", &TELRing::update)
      .def(
          "getblk",
          [](TELRing& x, int beam, int blk) {
            auto [buffer, size] = x.getblk(beam, blk);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "blk"_a)
      .def(
          "getblks",
          [](TELRing& x, int beam, int blk0, int blkN) {
            auto [buffer, size] = x.getblks(beam, blk0, blkN);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "blk0"_a, "blkN"_a)
      .def(
          "getslice",
          [](TELRing& x, int beam, double tbeg, double tend) {
            auto [buffer, size] = x.getslice(beam, tbeg, tend);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "tbeg"_a, "tend"_a);

  nb::class_<FRBRing>(m, "FRBRing")
      /** Constructor. **/
      .def(nb::init<>())

      /** Class properties. **/
      .def_prop_ro("mode",
                   [](FRBRing& x) {
                     switch (x.mode()) {
                       default:
                       case READ:
                         return "r";
                       case WRITE:
                         return "w";
                     }
                   })
      .def_prop_ro("opened", [](FRBRing& x) { return x.opened(); })

      /** PART I: Data properties. **/
      .def_prop_ro("nf", [](FRBRing& x) { return x.nf(); })
      .def_prop_ro("fh", [](FRBRing& x) { return x.fh(); })
      .def_prop_ro("fl", [](FRBRing& x) { return x.fl(); })
      .def_prop_ro("df", [](FRBRing& x) { return x.df(); })
      .def_prop_ro("bw", [](FRBRing& x) { return x.bw(); })
      .def_prop_ro("dt", [](FRBRing& x) { return x.dt(); })
      .def_prop_ro("nbits", [](FRBRing& x) { return x.nbits(); })
      .def_prop_ro("nstokes", [](FRBRing& x) { return x.nstokes(); })
      .def_prop_ro("flipped", [](FRBRing& x) { return x.flipped(); })

      /** PART II: Observation properties. **/
      .def_prop_ro("ra", [](FRBRing& x) { return x.ra(); })
      .def_prop_ro("dec", [](FRBRing& x) { return x.dec(); })
      .def_prop_ro("source", [](FRBRing& x) { return x.source(); })
      .def_prop_ro("beammode", [](FRBRing& x) { return x.beammode(); })
      .def_prop_ro("observer", [](FRBRing& x) { return x.observer(); })
      .def_prop_ro("gtaccode", [](FRBRing& x) { return x.gtaccode(); })
      .def_prop_ro("gtactitle", [](FRBRing& x) { return x.gtactitle(); })
      .def_prop_ro("antmaskpol1", [](FRBRing& x) { return x.antmaskpol1(); })
      .def_prop_ro("antmaskpol2", [](FRBRing& x) { return x.antmaskpol2(); })
      .def_prop_ro("antspol1", [](FRBRing& x) { return x.antspol1(); })
      .def_prop_ro("antspol2", [](FRBRing& x) { return x.antspol2(); })

      /** PART III: Beam steering and tiling properties. **/
      .def_prop_ro("beamid", [](FRBRing& x) { return x.beamid(); })
      .def_prop_ro("hostid", [](FRBRing& x) { return x.hostid(); })
      .def_prop_ro("nbeams", [](FRBRing& x) { return x.nbeams(); })
      .def_prop_ro("hostname", [](FRBRing& x) { return x.hostname(); })
      .def_prop_ro("npcbaselines", [](FRBRing& x) { return x.npcbaselines(); })
      .def_prop_ro("nbeamspernode", [](FRBRing& x) { return x.nbeamspernode(); })
      .def_prop_ro("beamras", [](FRBRing& x) { return x.beamras(); })
      .def_prop_ro("beamdecs", [](FRBRing& x) { return x.beamdecs(); })

      /** Summarise all properties as a dictionary. **/
      /** Summarise all properties as a dictionary. **/
      .def("header",
           [](FRBRing& x) {
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
      .def_prop_ro("maxblks", [](FRBRing& x) { return x.maxblks(); })
      .def_prop_ro("blksamps", [](FRBRing& x) { return x.blksamps(); })

      .def_prop_ro("empty", [](FRBRing& x) { return x.empty(); })
      .def_prop_ro("status", [](FRBRing& x) { return x.status(); })
      .def_prop_ro("active", [](FRBRing& x) { return x.active(); })

      .def_prop_ro("curblk", [](FRBRing& x) { return x.curblk(); })
      .def_prop_ro("currec", [](FRBRing& x) { return x.currec(); })
      .def_prop_ro("begblk", [](FRBRing& x) { return x.begblk(); })
      .def_prop_ro("endblk", [](FRBRing& x) { return x.endblk(); })

      .def_prop_ro("blksize", [](FRBRing& x) { return x.blksize(); })
      .def_prop_ro("size", [](FRBRing& x) { return x.size(); })

      .def_prop_ro("blktime", [](FRBRing& x) { return x.blktime(); })
      .def_prop_ro("curtime", [](FRBRing& x) { return x.curtime(); })
      .def_prop_ro("begtime", [](FRBRing& x) { return x.begtime(); })
      .def_prop_ro("endtime", [](FRBRing& x) { return x.endtime(); })
      .def("timeofblk", &FRBRing::timeofblk, "blk"_a)
      .def_prop_ro("timestamps", [](FRBRing& x) { return x.timestamps(); })

      /** Public methods. **/
      .def("open",
           [](FRBRing& x, std::string mode) {
             if (mode == "r") {
               x.open(READ);
             } else if ((mode == "w") || (mode == "rw")) {
               x.open(WRITE);
             } else {
               throw std::runtime_error("MODE DOESN'T EXIST. ABORT.");
             }
             return x;
           })
      .def("close", &FRBRing::close)
      .def("update", &FRBRing::update)
      .def(
          "getblk",
          [](FRBRing& x, int beam, int blk) {
            auto [buffer, size] = x.getblk(beam, blk);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "blk"_a)
      .def(
          "getblks",
          [](FRBRing& x, int beam, int blk0, int blkN) {
            auto [buffer, size] = x.getblks(beam, blk0, blkN);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "blk0"_a, "blkN"_a)
      .def(
          "getslice",
          [](FRBRing& x, int beam, double tbeg, double tend) {
            auto [buffer, size] = x.getslice(beam, tbeg, tend);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "tbeg"_a, "tend"_a)
      .def(
          "getburst",
          [](FRBRing& x, int beam, double t0, double dm, double width) {
            auto [buffer, size] = x.getburst(beam, t0, dm, width);
            size_t nf = x.nf();
            size_t nt = (size_t)(size / nf);
            return Array(buffer, {nt, nf},
                         nb::capsule(buffer, [](void* p) noexcept { delete[] (unsigned char*)p; }));
          },
          "beam"_a, "t0"_a, "dm"_a, "width"_a);
}
