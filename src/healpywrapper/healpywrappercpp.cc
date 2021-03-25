#include "ducc0/infra/communication.cc"
#include "ducc0/infra/string_utils.cc"
#include "ducc0/infra/system.cc"
#include "ducc0/infra/threading.cc"
#include "ducc0/infra/types.cc"
#include "ducc0/math/geom_utils.cc"
#include "ducc0/math/pointing.cc"
#include "ducc0/math/pointing.h"
#include "ducc0/math/space_filling.cc"
#include "ducc0/sharp/sharp.cc"
#include "ducc0/sharp/sharp_almhelpers.cc"
#include "ducc0/sharp/sharp_core.cc"
#include "ducc0/sharp/sharp_geomhelpers.cc"
#include "ducc0/sharp/sharp_ylmgen.cc"

#include <complex>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

#include "ducc0/bindings/pybind_utils.h"
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/string_utils.h"
#include "ducc0/math/constants.h"
#include "ducc0/sharp/sharp.h"
#include "ducc0/sharp/sharp_almhelpers.h"
#include "ducc0/sharp/sharp_geomhelpers.h"
#include "ducc0/healpix/healpix_base.h"
#include "ducc0/healpix/healpix_tables.h"
#include "ducc0/healpix/healpix_base.cc"
#include "ducc0/healpix/healpix_tables.cc"

using namespace ducc0;

namespace py = pybind11;

using a_d = py::array_t<double>;
using a_d_c = py::array_t<double, py::array::c_style | py::array::forcecast>;
using a_c_c =
    py::array_t<complex<double>, py::array::c_style | py::array::forcecast>;


/* Creates a geometry information describing a HEALPix map with an
   Nside parameter `nside`.
   The zbounds parameter `zbounds` is an array of two elements, with zbounds[0] <= zbounds[1].
   Every ring whose cos(theta) is not inside of the zbounds interval will not be included.
   `weight` contains the relative ring weights and must have 2*nside entries.
   If `weight` is a null pointer, all weights are assumed to be 1. */
unique_ptr<sharp_geom_info> sharp_make_zbounds_healpix_geom_info (size_t nside, ptrdiff_t stride, double* zbounds,
    const double *weight)
  {
  size_t npix=nside*nside*12;
  size_t nrings = 4*nside-1;
  size_t* rings = new size_t[nrings];

  size_t nrings_new = 0;

  for (size_t m = 0; m < nrings; ++m){
    size_t ring = m + 1;
    double cth;
    size_t northring = (ring>2*nside) ? 4*nside-ring : ring;
    if (northring < nside){
      cth = cos(2*asin(northring/(sqrt(6.)*nside)));
    } else {
      double fact1 = (8.*nside)/npix;
      cth = (2*nside-northring)*fact1;
    }
    if (northring != ring){ // southern hemisphere
      cth *= -1;
    }
    // if cth in zbounds, add ring to ring list
    if (cth >= zbounds[0] && cth <= zbounds[1]){
      rings[nrings_new] = ring;
      nrings_new ++;
    }
  }
   unique_ptr<sharp_geom_info> ginfo = ducc0::detail_sharp::sharp_make_subset_healpix_geom_info (nside, stride, nrings_new,
       rings, weight);
   delete rings;
   return ginfo;
  }

/* Returns the number of pixels of a map with parameters `nside` and `zbounds`.
  This number is lesser or equal to 12*nside**2. */
int get_npix(size_t nside, ptrdiff_t stride, a_d &zbounds){
  auto zb = zbounds.mutable_unchecked<1>();
  unique_ptr<sharp_geom_info> ginfo = sharp_make_zbounds_healpix_geom_info(nside, 1, &zb[0], NULL);
  int64_t npix = 0;
  for (long unsigned int i = 0; i < ginfo->nrings(); ++i){
    npix += ginfo->nph(i);
  }
  return npix;
}

/* Calculates the required offset of the input map for the 
  map2alm and alm2map routines, 
  when zbounds[1] is not equal to 1.
  It is equal to the number of pixels in the rings
  whose theta parameter is lesser than zbounds[1]. */
int offset(size_t nside, ptrdiff_t stride, a_d &zbounds){
  auto zb = zbounds.mutable_unchecked<1>();
  double * zb_ptr = &zb[0];
  size_t npix=nside*nside*12;
  size_t nrings = 4*nside-1;
  int offset = 0;
  int nph = 0; // number of pixels in current ring

  for (size_t m = 0; m < nrings; ++m){
    size_t ring = m + 1;
    double cth;
    size_t northring = (ring>2*nside) ? 4*nside-ring : ring;
    double fact1 = (8.*nside)/npix;
    if (northring < nside){
      cth = cos(2*asin(northring/(sqrt(6.)*nside)));
      nph = 4*northring;
    } else {
      cth = (2*nside-northring)*fact1;
      nph = 4*nside;
    }
    if (northring != ring){ // southern hemisphere
      cth *= -1;
    }
    if (cth >= zb_ptr[0] && cth <= zb_ptr[1]){
      return offset;
    } else {
      offset += nph;
    }
    }
  return 0;
  }

a_c_c map2alm(const a_d_c &map, const int64_t nside, const int64_t lmax,
                 const int64_t mmax, const int nthreads, a_d &zbounds) {

  // make triangular alm info
  MR_assert(mmax >= 0, "negative mmax");
  MR_assert(mmax <= lmax, "mmax must not be larger than lmax");
  unique_ptr<sharp_alm_info> ainfo =
      sharp_make_triangular_alm_info(lmax, mmax, 1);

  // set healpix zbounds geometry
  auto zb = zbounds.mutable_unchecked<1>();
  MR_assert(nside > 0, "bad nside value");
  unique_ptr<sharp_geom_info> ginfo = sharp_make_zbounds_healpix_geom_info(nside, 1, &zb[0], NULL);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto mr = map.unchecked<1>();
  auto ar = alm.mutable_unchecked<1>();

  sharp_map2alm(&ar[0], &mr[0], *ginfo, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}


void alm2map(const a_c_c &alm, const int64_t nside, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds, a_d_c &map){

    // make triangular alm info
    MR_assert(mmax >= 0, "negative mmax");
    MR_assert(mmax <= lmax, "mmax must not be larger than lmax");
    unique_ptr<sharp_alm_info> ainfo =
      sharp_make_triangular_alm_info(lmax, mmax, 1);

    // set healpix zbounds geometry
    auto zb = zbounds.mutable_unchecked<1>();
    MR_assert(nside > 0, "bad nside value");
    unique_ptr<sharp_geom_info> ginfo = sharp_make_zbounds_healpix_geom_info(nside, 1, &zb[0], NULL);

    int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
    MR_assert (alm.size()==n_alm,
        "incorrect size of a_lm array"); 

    auto mr=map.mutable_unchecked<1>();
    auto ar=alm.unchecked<1>();

    sharp_alm2map(&ar[0], &mr[0], *ginfo, *ainfo, 0, nthreads);
    return;
  }


/* binders */
using namespace pybind11;
PYBIND11_MODULE(healpywrappercpp, m) {
  auto my_submodule = m.def_submodule("sphtfunc");
  my_submodule.doc() =
      "Spherical Harmonic Transform functions module";
  my_submodule.def("map2alm", &map2alm, "map"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);
  my_submodule.def("alm2map", &alm2map, "map"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a, "map"_a);
  my_submodule.def("offset", &offset, "nside"_a, "stride"_a, "zbounds"_a);
  my_submodule.def("get_npix", &get_npix, "nside"_a, "stride"_a, "zbounds"_a);
}
