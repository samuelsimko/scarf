/*
 * Scarf binder file
 *
 *    Created 2021 by Samuel Simko.
 *    Based on the DUCC source code by Martin Reinecke.
 *
 * This file is subject to the terms and conditions of the MIT License.
 */

#include "docstrings.h"
#include "functions.h"

using namespace ducc0;
using namespace std;

namespace py = pybind11;

optional<long int> opt_int;
optional<a_d> opt_zb;

using namespace pybind11;

PYBIND11_MODULE(scarf, m) {

  m.attr("__version__") = "0.1.0";

  m.doc() = module_ds;

  /// SHT functions
  m.def("map2alm", &map2alm, map2alm_ds, "map"_a, "lmax"_a, "mmax"_a = opt_int,
        "nthreads"_a = 1, "zbounds"_a = opt_zb);
  m.def("map2alm_spin", &map2alm_spin, map2alm_ds, "map"_a, "spin"_a, "lmax"_a,
        "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb);
  m.def("alm2map", &alm2map, alm2map_ds, "alm"_a, "nside"_a, "lmax"_a,
        "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb);
  m.def("alm2map_spin", &alm2map_spin, alm2map_spin_ds, "alm"_a, "spin"_a,
        "nside"_a, "lmax"_a, "mmax"_a = opt_int, "nthreads"_a = 1,
        "zbounds"_a = opt_zb);
  m.def("GL_wg", &GL_wg, GL_wg_ds, "n"_a);
  m.def("GL_xg", &GL_xg, GL_xg_ds, "n"_a);

  /// Binders for the Geometry class
  py::class_<sharp_standard_geom_info>(m, "Geometry", geometry_ds)
      .def(py::init(&GeometryInformation), "nrings"_a, "nph"_a, "ofs"_a,
           "stride"_a, "phi0"_a, "theta"_a, "wgt"_a)

      /// attributes
      .def_readwrite("theta", &sharp_standard_geom_info::theta_array)
      .def_readwrite("phi0", &sharp_standard_geom_info::phi0_array)
      .def_readwrite("weight", &sharp_standard_geom_info::weight_array)
      .def_readwrite("nph", &sharp_standard_geom_info::nph_array)
      .def_readwrite("ofs", &sharp_standard_geom_info::ofs_array)

      /// methods 
      .def("get_nrings", &sharp_standard_geom_info::nrings, get_nrings_ds)
      .def("get_nph", &sharp_standard_geom_info::nph, get_nph_ds, "iring"_a)
      .def("get_nphmax", &sharp_standard_geom_info::nphmax, get_nphmax_ds)
      .def("get_theta", &sharp_standard_geom_info::theta, get_theta_ds,
           "iring"_a)
      .def("get_phi0", &sharp_standard_geom_info::phi0, get_phi0_ds, "iring"_a)
      .def("get_weight", &sharp_standard_geom_info::weight, get_weight_ds,
           "iring"_a)
      .def("get_ofs", &sharp_standard_geom_info::ofs, get_ofs_ds, "iring"_a)
      .def("map2alm", &map2alm_ginfo, map2alm_ginfo_ds, "map"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("alm2map", &alm2map_ginfo, alm2map_ginfo_ds, "alm"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("map2alm_spin", &map2alm_spin_ginfo, map2alm_spin_ginfo_ds, "map"_a,
           "spin"_a, "lmax"_a, "mmax"_a = opt_int, "nthreads"_a = 1,
           "zbounds"_a = opt_zb)
      .def("alm2map_spin", &alm2map_spin_ginfo, alm2map_spin_ds, "alm"_a,
           "spin"_a, "lmax"_a, "mmax"_a = opt_int, "nthreads"_a = 1,
           "zbounds"_a = opt_zb)
      .def("alm2phase", &alm2phase_ginfo, alm2phase_ds, "alm"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("phase2alm", &phase2alm_ginfo, phase2alm_ds, "phase"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("phase2map", &phase2map_ginfo, phase2map_ds, "phase"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("map2phase", &map2phase_ginfo, map2phase_ds, "map"_a, "lmax"_a,
           "mmax"_a = opt_int, "nthreads"_a = 1, "zbounds"_a = opt_zb)
      .def("alm2phase_spin", &alm2phase_spin_ginfo, alm2phase_spin_ds, "alm"_a,
           "spin"_a, "lmax"_a, "mmax"_a = opt_int, "nthreads"_a = 1,
           "zbounds"_a = opt_zb)
      .def("phase2alm_spin", &phase2alm_spin_ginfo, phase2alm_spin_ds,
           "phase"_a, "spin"_a, "lmax"_a, "mmax"_a = opt_int, "nthreads"_a = 1,
           "zbounds"_a = opt_zb);

  /// Geometry functions
  m.def("healpix_geometry", &sharp_make_healpix_geom_info, healpix_geometry_ds,
        "nside"_a, "stride"_a);

  m.def("gauss_geometry", &gauss_geometry, gauss_geometry_ds, "nside"_a,
        "stride"_a);
}
