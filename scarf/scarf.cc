#include "ducc0/infra/communication.cc"
#include "ducc0/infra/system.cc"
#include "ducc0/infra/threading.cc"
#include "ducc0/infra/types.cc"
#include "ducc0/math/geom_utils.cc"
#include "ducc0/math/pointing.cc"
#include "ducc0/math/pointing.h"
#include "ducc0/math/space_filling.cc"

#include <complex>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>
#include "ducc0/bindings/pybind_utils.h"
#include "ducc0/healpix/healpix_base.h"
#include "ducc0/healpix/healpix_tables.h"
#include "ducc0/healpix/healpix_base.cc"
#include "ducc0/healpix/healpix_tables.cc"

#include "ducc0/sharp/sht.h"
#include "ducc0/sharp/sht.cc"
#include "ducc0/sharp/sharp.h"
#include "ducc0/sharp/sharp.cc"
#include "ducc0/sharp/sharp_geomhelpers.h"
#include "ducc0/sharp/sharp_geomhelpers.cc"
#include "ducc0/sharp/sharp_almhelpers.h"
#include "ducc0/sharp/sharp_almhelpers.cc"
#include "ducc0/infra/string_utils.h"
#include "ducc0/infra/string_utils.cc"
#include "ducc0/infra/error_handling.h"
#include "ducc0/math/constants.h"

// #include "phase.h"

using namespace ducc0;
using namespace std;

namespace py = pybind11;

using a_d = py::array_t<double>;
using a_s = py::array_t<size_t>;
using a_li = py::array_t<long int>;
using a_d_c = py::array_t<double, py::array::c_style | py::array::forcecast>;
using a_c_c =
    py::array_t<complex<double>, py::array::c_style | py::array::forcecast>;

unique_ptr<sharp_alm_info> set_triangular_alm_info (int64_t lmax, int64_t mmax)
  {
  MR_assert(mmax>=0,"negative mmax");
  MR_assert(mmax<=lmax,"mmax must not be larger than lmax");
  return sharp_make_triangular_alm_info(lmax,mmax,1);
  }


using namespace ducc0::detail_sharp;

sharp_standard_geom_info* GeometryInformation(size_t nrings, a_s &nph, a_li &ofs, ptrdiff_t stride, a_d &phi0, a_d &theta, a_d &wgt){
  auto nph_p = nph.unchecked<1>();
  auto ofs_p = ofs.unchecked<1>();
  auto phi0_p = phi0.unchecked<1>();
  auto theta_p = theta.unchecked<1>();
  auto wgt_p = wgt.unchecked<1>();
  auto temp_p = make_unique<sharp_standard_geom_info>(nrings, &nph_p[0], &ofs_p[0], stride, &phi0_p[0], &theta_p[0], &wgt_p[0]);
  return temp_p.release();
}

// creates a new geometry info, keeping the rings only in zbounds
sharp_standard_geom_info * keep_rings_in_zbounds(sharp_standard_geom_info &ginfo, double * zbounds){
   size_t nrings = ginfo.nrings();

   vector<size_t> nph;
   vector<ptrdiff_t> ofs;
   ptrdiff_t stride = 1;
   vector<double> phi0;
   vector<double> theta;
   vector<double> wgt;

  size_t iring = 0;
  size_t nrings_new = 0;
  for (; iring < nrings; ++iring){
    if (cos(ginfo.theta(iring)) >= zbounds[0] && cos(ginfo.theta(iring)) <= zbounds[1]){
        nrings_new += 1;
        // add ring info to new structure
        nph.push_back(ginfo.nph(iring));
        ofs.push_back(ginfo.ofs(iring));
        phi0.push_back(ginfo.phi0(iring));
        theta.push_back(ginfo.theta(iring));
        wgt.push_back(ginfo.weight(iring));
        }
  }
  return make_unique<sharp_standard_geom_info>(nrings_new, &nph[0], &ofs[0], stride, &phi0[0], &theta[0], &wgt[0]).release();
}

sharp_standard_geom_info * sharp_make_standard_healpix_geom_info(size_t nside, size_t stride){
  return sharp_make_healpix_geom_info(nside, stride).release();
}


a_c_c map2alm_ginfo(sharp_standard_geom_info *ginfo, a_d_c map, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_standard_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

  // make triangular alm info
  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto mr = map.unchecked<1>();
  auto ar = alm.mutable_unchecked<1>();

  sharp_map2alm(&ar[0], &mr[0], *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

a_c_c map2alm_spin_ginfo(sharp_standard_geom_info *ginfo, const a_d_c &map, int64_t spin,
    const int64_t lmax, const int64_t mmax, const int nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_standard_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

  // make triangular alm info
  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(vector<size_t>{2,size_t(n_alm)});
  auto ar = alm.mutable_unchecked<2>();

  int64_t npix = 0;
  for (size_t i = 0; i < ginfo->nrings(); ++i){;
    npix += ginfo->nph(i);
  }

  auto mr=map.unchecked<2>();
  MR_assert ((mr.shape(0)==2)&&(mr.shape(1)==npix), "incorrect size of map array");

  sharp_map2alm_spin(spin, &ar(0,0), &ar(1,0), &mr(0, 0), &mr(1, 0), *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

a_d_c alm2map_ginfo(sharp_standard_geom_info *ginfo, const a_c_c &alm, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

    auto zb = zbounds.mutable_unchecked<1>();
    sharp_standard_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

    // make triangular alm info
    unique_ptr<sharp_alm_info> ainfo =
      set_triangular_alm_info (lmax, mmax);

    int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
    MR_assert (alm.size()==n_alm, "incorrect size of a_lm array"); 

    size_t npix = 0;
    for (size_t i = 0; i < ginfo->nrings(); ++i){;
      npix += ginfo->nph(i);
    }

    a_d_c map(npix, 0);
    auto mr=map.mutable_unchecked<1>();

    for (size_t i = 0; i < npix; ++i){
      mr[i] = 0;
    }

    auto ar=alm.unchecked<1>();

    sharp_alm2map(&ar[0], &mr[0], *ginfo_new, *ainfo, 0, nthreads);
    return map;
}

a_d_c alm2map_spin_ginfo(sharp_standard_geom_info *ginfo, const a_c_c &alm, int64_t spin,
    const int64_t lmax, const int64_t mmax, const int nthreads, a_d &zbounds){

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_standard_geom_info* ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

  // make triangular alm info
  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  auto ar=alm.unchecked<2>();
  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  MR_assert((ar.shape(0)==2)&&(ar.shape(1)==n_alm),
    "incorrect size of a_lm array");

  int64_t npix = 0;
  for (size_t i = 0; i < ginfo->nrings(); ++i){;
    npix += ginfo->nph(i);
  }

  a_d_c map(vector<size_t>{2,size_t(npix)});
  auto mr=map.mutable_unchecked<2>();
  sharp_alm2map_spin(spin, &ar(0,0), &ar(1,0), &mr(0, 0), &mr(1, 0), *ginfo_new, *ainfo, 0, nthreads);
  return map;
  }

a_c_c map2alm(const a_d_c &map, const int64_t lmax,
                 const int64_t mmax, const int nthreads, a_d &zbounds) {

    const int64_t npix = map.ndim() == 1 ? map.size() : map.shape(1);
    const int64_t nside = (int) sqrt(npix/12);
    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_ginfo(ginfo, map, lmax, mmax, nthreads, zbounds);

}

a_c_c map2alm_spin(const a_d_c &map, int64_t spin, const int64_t lmax,
                 const int64_t mmax, const int nthreads, a_d &zbounds) {

    const int64_t npix = map.ndim() == 1 ? map.size() : map.shape(1);
    const int64_t nside = (int) sqrt(npix/12);
    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_spin_ginfo(ginfo, map, spin, lmax, mmax, nthreads, zbounds);
}

a_d_c alm2map(const a_c_c &alm, const int64_t nside, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return alm2map_ginfo(ginfo, alm, lmax, mmax, nthreads, zbounds);
  }


a_d_c alm2map_spin(const a_c_c &alm, int64_t spin, const int64_t nside, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

  sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
  return alm2map_spin_ginfo(ginfo, alm, spin, lmax, mmax, nthreads, zbounds);
  }

py::array GL_wg(size_t n)
  {
  auto res = make_Pyarr<double>({n});
  auto res2 = to_mav<double,1>(res, true);
  GL_Integrator integ(n);
  auto wgt = integ.weights();
  for (size_t i=0; i<res2.shape(0); ++i)
    res2.v(i) = wgt[i];
  return move(res);
  }

py::array GL_xg(size_t n)
  {
  auto res = make_Pyarr<double>({n});
  auto res2 = to_mav<double,1>(res, true);
  GL_Integrator integ(n);
  auto x = integ.coords();
  for (size_t i=0; i<res2.shape(0); ++i)
    res2.v(i) = -x[i];
  return move(res);
  }


/* binders */

using namespace pybind11;
PYBIND11_MODULE(scarf, m) {
  m.attr("__version__") = "0.1.0";

  m.doc() = R"pbdoc(
  Spherical harmonics transform library for CMB lensing
  )pbdoc";

  m.def("map2alm", &map2alm, R"pbdoc(
  Computes alms from a temperature map.

  Parameters
  ----------
  map : np.array, shape (:math:`N_{pix}`,)
    The temperature map
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{alm}`)
    Temperature alm
  )pbdoc", "map"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a);


  m.def("map2alm_spin", &map2alm_spin, R"pbdoc(
  Computes alms from a polarisation map.

  Parameters
  ----------
  map : np.array, shape (:math:`N_{pix}`,)
    The polarisation map.
  spin : int, scalar
    The spin of the field.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{alm}`)
    Polarisation alm
  )pbdoc", "map"_a, "spin"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("alm2map", &alm2map, R"pbdoc(
  Computes a Healpix temperature map from alm.

  Parameters
  ----------
  alm : np.array, shape (:math:`N_{alm}`,)
    The temperature alm.
  nside : int, scalar
    The number of sides for the output map. This is the number of pixels along the diagonal of a healpix-grid-base-pixel.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{pix}`)
    Temperature map
  )pbdoc", "alm"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("alm2map_spin", &alm2map_spin, R"pbdoc(
  Computes a Healpix polarisation map from alm.

  Parameters
  ----------
  alm : np.array, shape (:math:`N_{alm}`,)
    The polarisation alm
  spin : int, scalar
    The spin of the field.
  nside : int, scalar
    The number of sides for the output map. This is the number of pixels along the diagonal of a healpix-grid-base-pixel.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{pix}`)
    Temperature map
  )pbdoc", "alm"_a, "spin"_a, "nside"_a, 
      "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("GL_wg", &GL_wg,  R"pbdoc(
  Computes Gauss-Legendre quadrature weights.

  Parameters
  ----------
  n: int, scalar
    The number of GL sample points used (integrates exactly degree :math:`2n - 1` polynomials; typically :math:`n=\ell_{\rm max} + 1`)

  Returns
  -------
  np.array, shape (:math:`n`)
  )pbdoc","n"_a);

  m.def("GL_xg", &GL_xg,  R"pbdoc(
  Computes Gauss-Legendre quadrature sample points. Output is ordered from 1 to -1.

  Parameters
  ----------
  n: int, scalar
    The number of GL sample points used (integrates exactly degree :math:`2n - 1` polynomials; typically :math:`n=\ell_{\rm max} + 1`)

  Returns
  -------
  np.array, shape (:math:`n`)
  )pbdoc","n"_a);

  py::class_<sharp_standard_geom_info>(m ,"Geometry", R"pbdoc(
  Creates a geometry specified by the user.

  Parameters
  ----------
  nrings : int, scalar
    The number of rings
  nph : int,  shape (:math:`N_{rings}`)
    TBD
  ofs : int,  shape (:math:`N_{rings}`)
    TBD
  stride : int, scalar
    The stride between two consecutive pixels on the map
  phi0 : float, shape (:math:`N_{rings}`)
    TBD
  theta : float,  shape (:math:`N_{rings}`)
    The latitude angle of each ring
  wgt : float,  shape (:math:`N_{rings}`)
    The weighting for each ring
  
  Returns
  -------
  Geometry
    A geometry as specified by the user
  )pbdoc")
    .def(py::init(&GeometryInformation), "nrings"_a, "nph"_a, "ofs"_a, "stride"_a, "phi0"_a, "theta"_a, "wgt"_a )
    .def("get_nrings", &sharp_standard_geom_info::nrings, R"pbdoc(
    Returns the number of rings of the geometry.

    Parameters
    ----------

    Returns
    -------
    int
      The number of rings of the geometry
      
    )pbdoc")
    .def("get_nph", &sharp_standard_geom_info::nph, R"pbdoc(
    Returns the number of pixels in the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    int
      The number of pixels in the ring.
      
    )pbdoc", "iring"_a)
    .def("get_nphmax", &sharp_standard_geom_info::nphmax,  R"pbdoc(
    Returns the maximum number of pixels of a ring.

    Parameters
    ----------

    Returns
    -------
    int
      The maximum number of pixels of a ring.
      
    )pbdoc")
    .def("get_theta", &sharp_standard_geom_info::theta,  R"pbdoc(
    Returns the lattitude, in radiants, of a specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The latitude of the ring
      
    )pbdoc", "iring"_a)
    .def("get_cth", &sharp_standard_geom_info::cth,  R"pbdoc(
    Returns the cosinus of the latitude of a specified ring. 

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The cosinus of the latitude of the ring.
      
    )pbdoc", "iring"_a)
    .def("get_sth", &sharp_standard_geom_info::sth,  R"pbdoc(
    Returns the sinus of the latitude of a specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The sinus of the latitude of the ring.
      
    )pbdoc", "iring"_a)
    .def("get_phi0", &sharp_standard_geom_info::phi0,  R"pbdoc(
    Returns the longitude of the first pixel of a ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The longitude of the first pixel of the ring.
      
    )pbdoc", "iring"_a)
    .def("get_weight", &sharp_standard_geom_info::weight,  R"pbdoc(
    Returns the weight of the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The weight of the ring. 
      
    )pbdoc", "iring"_a)
    .def("get_ofs", &sharp_standard_geom_info::ofs,  R"pbdoc(
    Returns the pixel offset of the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The pixel offset of the ring. 
      
    )pbdoc", "iring"_a)
    .def_readwrite("theta", &sharp_standard_geom_info::theta_array)
    .def_readwrite("phi0", &sharp_standard_geom_info::phi0_array)
    .def_readwrite("weight", &sharp_standard_geom_info::weight_array)
    .def_readwrite("cth", &sharp_standard_geom_info::cth_array)
    .def_readwrite("sth", &sharp_standard_geom_info::sth_array)
    .def_readwrite("nph", &sharp_standard_geom_info::nph_array)
    .def_readwrite("ofs", &sharp_standard_geom_info::ofs_array)
    //.def("pair", &sharp_standard_geom_info::pair, "iring"_a)
    //.def("get_ring", &sharp_standard_geom_info::get_ring, "weighted"_a, "iring"_a, "map"_a, "ringtmp"_a)
    //.def("add_ring", &sharp_standard_geom_info::add_ring, "weighted"_a, "iring"_a, "ringtmp"_a, "map"_a)
    .def("map2alm", &map2alm_ginfo, "map"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("alm2map", &alm2map_ginfo, "alm"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("map2alm_spin", &map2alm_spin_ginfo, "map"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("alm2map_spin", &alm2map_spin_ginfo, "alm"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("healpix_geometry", &sharp_make_standard_healpix_geom_info, R"pbdoc(
  Creates a HEALPix geometry given the nside.

  Parameters
  ----------
  nside : int, scalar
    The nside of the HEALPix geometry
  stride: int, scalar
    The stride between two consecutive pixels on the map

  Returns
  -------
  Geometry
    A Scarf geometry following the HEALPix scheme
  )pbdoc", "nside"_a, "stride"_a);

}
