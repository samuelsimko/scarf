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
#include <pybind11/stl.h>
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

#include "phase.h"

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
sharp_standard_geom_info * keep_rings_in_zbounds(sharp_standard_geom_info &ginfo, a_d &zbounds){

  auto zb = zbounds.mutable_unchecked<1>();

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
    if (cos(ginfo.theta(iring)) >= zb[0] && cos(ginfo.theta(iring)) <= zb[1]){
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

a_c_c map2alm_ginfo(sharp_standard_geom_info *ginfo, a_d_c map, size_t lmax, optional<size_t> mmax_opt, size_t nthreads, a_d &zbounds) {

  size_t mmax =  (mmax_opt.has_value()) ? mmax_opt.value() : lmax;

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

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
    const int64_t lmax, optional<int64_t> mmax_opt, const int nthreads, a_d &zbounds) {


  size_t mmax =  (mmax_opt.has_value()) ? mmax_opt.value() : lmax;

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

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
      optional<int64_t> mmax_opt, const int nthreads, a_d &zbounds){

  size_t mmax =  (mmax_opt.has_value()) ? mmax_opt.value() : lmax;

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

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
    const int64_t lmax, optional<int64_t> mmax_opt, const int nthreads, a_d &zbounds){

  size_t mmax =  (mmax_opt.has_value()) ? mmax_opt.value() : lmax;
  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

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
    const optional<int64_t> mmax, const int nthreads, a_d &zbounds) {

    const int64_t npix = map.ndim() == 1 ? map.size() : map.shape(1);
    const int64_t nside = (int) sqrt(npix/12);
    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_ginfo(ginfo, map, lmax, mmax, nthreads, zbounds);

}

a_c_c map2alm_spin(const a_d_c &map, int64_t spin, const int64_t lmax,
    const optional<int64_t> mmax, const int nthreads, a_d &zbounds) {

    const int64_t npix = map.ndim() == 1 ? map.size() : map.shape(1);
    const int64_t nside = (int) sqrt(npix/12);
    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_spin_ginfo(ginfo, map, spin, lmax, mmax, nthreads, zbounds);
}

a_d_c alm2map(const a_c_c &alm, const int64_t nside, const int64_t lmax,
    const optional<int64_t> mmax, const int nthreads, a_d &zbounds){

    sharp_standard_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return alm2map_ginfo(ginfo, alm, lmax, mmax, nthreads, zbounds);
  }


a_d_c alm2map_spin(const a_c_c &alm, int64_t spin, const int64_t nside, const int64_t lmax,
    const optional<int64_t> mmax, const int nthreads, a_d &zbounds){

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

/* phase functions */

a_c_c alm2phase_ginfo(sharp_standard_geom_info *ginfo, const a_c_c &alm, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds, py::object &out){

    sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

    // make triangular alm info
    unique_ptr<sharp_alm_info> ainfo =
      set_triangular_alm_info (lmax, mmax);

    int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
    MR_assert (alm.size()==n_alm, "incorrect size of a_lm array"); 

    size_t npix = 0;
    for (size_t i = 0; i < ginfo->nrings(); ++i){;
      npix += ginfo->nph(i);
    }

    long unsigned int nchunks;
    long unsigned int chunksize;
    get_singular_chunk_info(ginfo_new->npairs(), (0==0) ? 128 : 64, 
        nchunks,chunksize);

    // auto phase_a2p = get_optional_Pyarr<complex<double>>(out, {2*chunksize, unsigned(mmax)+1, 1});
    auto phase_a2p = get_optional_Pyarr<complex<double>>(out, {1, 2*chunksize, unsigned(mmax)+1});
    auto phase_mav_a2p = to_mav<complex<double>, 3>(phase_a2p, true);

    py::buffer_info buf = phase_a2p.request();
    auto *ptr = static_cast<complex<double> *>(buf.ptr);

    // for safety reasons 
    for (size_t i = 0; i < 2*chunksize*(mmax+1)*1; ++i){
      ptr[i] = std::complex<double>(0);
    }

    auto ar=alm.unchecked<1>();

    sharp_alm2phase(&ar[0], phase_mav_a2p, *ginfo_new, *ainfo, 0, nthreads);
    return phase_a2p;
}


a_c_c phase2alm_ginfo(sharp_standard_geom_info *ginfo, a_c_c &phase_p2a, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto ar = alm.mutable_unchecked<1>();
  auto phase_mav_p2a = to_mav<complex<double>, 3>(phase_p2a, true); // second param false or true ?

  sharp_phase2alm(&ar[0], phase_mav_p2a, *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

a_d_c phase2map_ginfo(sharp_standard_geom_info *ginfo, a_c_c &phase, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  size_t npix = 0;
  for (size_t i = 0; i < ginfo->nrings(); ++i){;
    npix += ginfo->nph(i);
  }

  if (phase.shape(2) == 2){
    int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
    a_c_c alm(vector<size_t>{2, size_t(n_alm)});
    auto ar = alm.mutable_unchecked<2>();

    auto phase_mav = to_mav<complex<double>, 3>(phase, true);
    a_d_c map(vector<size_t>{2,size_t(npix)});
    auto mr=map.mutable_unchecked<2>();

    phase_job job(SHARP_Y, 2, {&ar(0, 0), &ar(1, 0)}, {&mr(0, 0), &mr(1, 0)}, phase_mav, *ginfo_new, *ainfo, 0, nthreads);
    phase_execute_phase2map(job, phase_mav, *ginfo_new, mmax, 2);
    return map;
  }

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto ar = alm.mutable_unchecked<1>();
  auto phase_mav = to_mav<complex<double>, 3>(phase, true);

  a_d_c map(npix);
  auto mr=map.mutable_unchecked<1>();

  for (size_t i = 0; i < npix; ++i){
    mr[i] = 0;
  }

  phase_job job(SHARP_Y, 0, {&ar[0]}, {&mr[0]}, phase_mav, *ginfo_new, *ainfo, 0, nthreads);
  phase_execute_phase2map(job, phase_mav, *ginfo_new, mmax, 0);
  return map;
}

a_c_c map2phase_ginfo(sharp_standard_geom_info *ginfo, a_d_c &map, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds,
    py::object &out) {

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  auto mr=map.mutable_unchecked<1>();

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto ar = alm.mutable_unchecked<1>();


  long unsigned int nchunks;
  long unsigned int chunksize;
  get_singular_chunk_info(ginfo_new->npairs(), (0==0) ? 128 : 64, 
      nchunks,chunksize);

  auto phase_m2p = get_optional_Pyarr<complex<double>>(out, {unsigned(map.ndim()), 2*chunksize, unsigned(mmax)+1});
  auto phase_mav = to_mav<complex<double>, 3>(phase_m2p, true);

  py::buffer_info buf = phase_m2p.request();
  auto *ptr = static_cast<complex<double> *>(buf.ptr);

  // for safety reasons
  for (size_t i = 0; i < 2*chunksize*(mmax+1)*1; ++i){
    ptr[i] = std::complex<double>(0);
  }

  phase_job job(SHARP_Yt, 0, {&ar[0]}, {&mr[0]}, phase_mav, *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  phase_execute_map2phase(job, phase_mav, *ginfo_new, mmax, 0);
  return phase_m2p;
}

a_c_c alm2phase_spin_ginfo(sharp_standard_geom_info *ginfo, const a_c_c &alm, const size_t spin, const int64_t lmax,
     const int64_t mmax, const int nthreads, a_d &zbounds, py::object &out){

    sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

    // make triangular alm info
    unique_ptr<sharp_alm_info> ainfo =
      set_triangular_alm_info (lmax, mmax);

    int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);

    size_t npix = 0;
    for (size_t i = 0; i < ginfo->nrings(); ++i){;
      npix += ginfo->nph(i);
    }

    long unsigned int nchunks;
    long unsigned int chunksize;
    get_singular_chunk_info(ginfo_new->npairs(), (0==0) ? 128 : 64, 
        nchunks,chunksize);

    auto phase = get_optional_Pyarr<complex<double>>(out, {2*chunksize, unsigned(mmax)+1, 2});
    auto phase_mav = to_mav<complex<double>, 3>(phase, true);

    auto ar=alm.unchecked<2>();
    MR_assert((ar.shape(0)==2)&&(ar.shape(1)==n_alm),
      "incorrect size of a_lm array");

    sharp_alm2phase_spin(spin, &ar(0, 0), &ar(1, 0), phase_mav, *ginfo_new, *ainfo, 0, nthreads);
    return phase;
}

a_c_c phase2alm_spin_ginfo(sharp_standard_geom_info *ginfo, a_c_c &phase, size_t spin, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  sharp_standard_geom_info *ginfo_new = (zbounds.is_none()) ? ginfo : keep_rings_in_zbounds(*ginfo, zbounds);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(vector<size_t>{2,size_t(n_alm)});
  auto ar = alm.mutable_unchecked<2>();
  auto phase_mav = to_mav<complex<double>, 3>(phase);

  sharp_phase2alm_spin(spin, &ar(0, 0), &ar(1, 0), phase_mav, *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

sharp_standard_geom_info * gauss_geometry(int64_t nrings, int64_t nphi)
      {
      MR_assert((nrings>0)&&(nphi>0),"bad grid dimensions");
      sharp_standard_geom_info * ginfo = sharp_make_2d_geom_info (nrings, nphi, 0., 1, nphi, "GL").release();
      return ginfo;
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
  )pbdoc", "map"_a, "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none());


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
                   "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none());

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
                   "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none());

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
          "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none());

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
    int
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
    .def("map2alm", &map2alm_ginfo, "map"_a, "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none())
    .def("alm2map", &alm2map_ginfo, "alm"_a, "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none())
    .def("map2alm_spin", &map2alm_spin_ginfo, "map"_a, "spin"_a, "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none())
    .def("alm2map_spin", &alm2map_spin_ginfo, "alm"_a, "spin"_a, "lmax"_a, "mmax"_a=py::none(), "nthreads"_a=1, "zbounds"_a=py::none())
    .def("alm2phase", &alm2phase_ginfo, "alm"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a, "out"_a=py::none())
    .def("phase2alm", &phase2alm_ginfo, "phase"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("phase2map", &phase2map_ginfo, "phase"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("map2phase", &map2phase_ginfo, "map"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a, "out"_a=py::none())
    .def("alm2phase_spin", &alm2phase_spin_ginfo, "alm"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a, "out"_a=py::none())
    .def("phase2alm_spin", &phase2alm_spin_ginfo, "phase"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a);


  m.def("healpix_geometry", &sharp_make_healpix_geom_info, R"pbdoc(
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

  m.def("gauss_geometry", &gauss_geometry, R"pbdoc(
  Creates a gauss geometry given nrings and nphi.

  Parameters
  ----------
  nrings : int, scalar
    The number of rings of the geometry
  nphi: int, scalar
    The number of pixels in each ring

  Returns
  -------
  Geometry
    A Scarf geometry following the gauss scheme
  )pbdoc", "nside"_a, "stride"_a);
}
