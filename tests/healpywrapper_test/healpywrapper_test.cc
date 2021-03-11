#include "ducc0/healpix/healpix_base.cc"
#include "ducc0/healpix/healpix_tables.cc"
#include "ducc0/infra/communication.cc"
#include "ducc0/infra/string_utils.cc"
#include "ducc0/infra/system.cc"
#include "ducc0/infra/threading.cc"
#include "ducc0/infra/types.cc"
#include "ducc0/math/geom_utils.cc"
#include "ducc0/math/pointing.cc"
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

using namespace ducc0;

namespace py = pybind11;

using a_d = py::array_t<double>;
using a_d_c = py::array_t<double, py::array::c_style | py::array::forcecast>;
using a_c_c =
    py::array_t<complex<double>, py::array::c_style | py::array::forcecast>;

a_c_c my_map2alm(const a_d_c &map, const int64_t nside, const int64_t lmax,
                 const int64_t mmax, const int nthreads) {
  // make triangular alm info
  MR_assert(mmax >= 0, "negative mmax");
  MR_assert(mmax <= lmax, "mmax must not be larger than lmax");
  unique_ptr<sharp_alm_info> ainfo =
      sharp_make_triangular_alm_info(lmax, mmax, 1);

  // set healpix geometry
  MR_assert(nside > 0, "bad nside value");
  int64_t npix = 12 * nside * nside;
  unique_ptr<sharp_geom_info> ginfo = sharp_make_healpix_geom_info(nside, 1);

  // calculate alms
  MR_assert(map.size() == npix, "incorrect size of map array");
  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto mr = map.unchecked<1>();
  auto ar = alm.mutable_unchecked<1>();

  sharp_map2alm(&ar[0], &mr[0], *ginfo, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

using namespace pybind11;
PYBIND11_MODULE(healpywrapper_test, m) {
  auto my_submodule = m.def_submodule("my_submodule");
  my_submodule.doc() =
      "My submodule, where I'll be testing the functions I implement.";
  my_submodule.def("my_map2alm", &my_map2alm, "map"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a);
}
