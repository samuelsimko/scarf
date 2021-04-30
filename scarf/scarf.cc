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
      // cth = 1 - (northring*northring)/(3*nside*nside); 
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
      // cth = 1 - (northring*northring)/(3*nside*nside); 
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



using namespace ducc0::detail_sharp;

sharp_geom_info* GeometryInformation(size_t nrings, a_s &nph, a_li &ofs, ptrdiff_t stride, a_d &phi0, a_d &theta, a_d &wgt){
  auto nph_p = nph.unchecked<1>();
  auto ofs_p = ofs.unchecked<1>();
  auto phi0_p = phi0.unchecked<1>();
  auto theta_p = theta.unchecked<1>();
  auto wgt_p = wgt.unchecked<1>();
  auto temp_p = make_unique<sharp_standard_geom_info>(nrings, &nph_p[0], &ofs_p[0], stride, &phi0_p[0], &theta_p[0], &wgt_p[0]);
  return temp_p.release();
}

// creates a new geometry info, keeping the rings only in zbounds
sharp_geom_info * keep_rings_in_zbounds(sharp_geom_info &ginfo, double * zbounds){
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



a_c_c map2alm_ginfo(sharp_geom_info *ginfo, a_d_c map, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_d_c map2alm_spin_ginfo(sharp_geom_info *ginfo, const a_d_c &map, int64_t spin,
    const int64_t lmax, const int64_t mmax, const int nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_d_c alm2map_ginfo(sharp_geom_info *ginfo, const a_c_c &alm, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

    auto zb = zbounds.mutable_unchecked<1>();
    sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_d_c alm2map_spin_ginfo(sharp_geom_info *ginfo, const a_c_c &alm, int64_t spin,
    const int64_t lmax, const int64_t mmax, const int nthreads, a_d &zbounds){

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info* ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_c_c map2alm(const a_d_c &map, const int64_t nside, const int64_t lmax,
                 const int64_t mmax, const int nthreads, a_d &zbounds) {

    sharp_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_ginfo(ginfo, map, lmax, mmax, nthreads, zbounds);

}

a_d_c map2alm_spin(const a_d_c &map, int64_t spin, const int64_t nside, const int64_t lmax,
                 const int64_t mmax, const int nthreads, a_d &zbounds) {

    sharp_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return map2alm_spin_ginfo(ginfo, map, spin, lmax, mmax, nthreads, zbounds);
}

a_d_c alm2map(const a_c_c &alm, const int64_t nside, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

    sharp_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
    return alm2map_ginfo(ginfo, alm, lmax, mmax, nthreads, zbounds);
  }


a_d_c alm2map_spin(const a_c_c &alm, int64_t spin, const int64_t nside, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds){

  sharp_geom_info *ginfo = sharp_make_healpix_geom_info(nside, 1).release();
  return alm2map_spin_ginfo(ginfo, alm, spin, lmax, mmax, nthreads, zbounds);
  }

/* phase functions */

a_c_c alm2phase_ginfo(sharp_geom_info *ginfo, const a_c_c &alm, const int64_t lmax,
      const int64_t mmax, const int nthreads, a_d &zbounds, py::object &out){

    auto zb = zbounds.mutable_unchecked<1>();
    sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

    auto phase_a2p = get_optional_Pyarr<complex<double>>(out, {2*chunksize, unsigned(mmax)+1, 1});
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

a_c_c phase2alm_ginfo(sharp_geom_info *ginfo, a_c_c &phase_p2a, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(n_alm);
  auto ar = alm.mutable_unchecked<1>();
  auto phase_mav_p2a = to_mav<complex<double>, 3>(phase_p2a); // second param false or true ?

  sharp_phase2alm(&ar[0], phase_mav_p2a, *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

a_d_c phase2map_ginfo(sharp_geom_info *ginfo, a_c_c &phase, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_c_c map2phase_ginfo(sharp_geom_info *ginfo, a_d_c &map, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds,
    py::object &out) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

  auto phase_m2p = get_optional_Pyarr<complex<double>>(out, {2*chunksize, unsigned(mmax)+1, unsigned(map.ndim())});
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

a_c_c alm2phase_spin_ginfo(sharp_geom_info *ginfo, const a_c_c &alm, const size_t spin, const int64_t lmax,
     const int64_t mmax, const int nthreads, a_d &zbounds, py::object &out){

    auto zb = zbounds.mutable_unchecked<1>();
    sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

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

a_c_c phase2alm_spin_ginfo(sharp_geom_info *ginfo, a_c_c &phase, size_t spin, size_t lmax, size_t mmax, size_t nthreads, a_d &zbounds) {

  auto zb = zbounds.mutable_unchecked<1>();
  sharp_geom_info *ginfo_new = keep_rings_in_zbounds(*ginfo, &zb[0]);

  unique_ptr<sharp_alm_info> ainfo =
    set_triangular_alm_info (lmax, mmax);

  int64_t n_alm = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
  a_c_c alm(vector<size_t>{2,size_t(n_alm)});
  auto ar = alm.mutable_unchecked<2>();
  auto phase_mav = to_mav<complex<double>, 3>(phase);

  sharp_phase2alm_spin(spin, &ar(0, 0), &ar(1, 0), phase_mav, *ginfo_new, *ainfo, SHARP_USE_WEIGHTS, nthreads);
  return alm;
}

/* binders */

using namespace pybind11;
PYBIND11_MODULE(scarf, m) {
  m.attr("__version__") = "0.1.0";

  m.doc() = R"pbdoc(
  Spherical harmonics transform library for CMB lensing
  )pbdoc";

  m.def("map2alm", &map2alm, R"pbdoc(
  Computes alms from a given map
  )pbdoc", "map"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);


  m.def("map2alm_spin", &map2alm_spin, R"pbdoc(
  Computes alms from a given map
  )pbdoc", "map"_a, "spin"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("alm2map", &alm2map, R"pbdoc(
    Computes a Healpix map given the alm.
  )pbdoc", "alm"_a, "nside"_a, "lmax"_a,
                   "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("alm2map_spin", &alm2map_spin, "alm"_a, "spin"_a, "nside"_a,
      "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a);

  m.def("offset", &offset, "nside"_a, "stride"_a, "zbounds"_a);
  m.def("get_npix", &get_npix, "nside"_a, "stride"_a, "zbounds"_a);

  py::class_<sharp_geom_info>(m ,"Geometry")
    .def(py::init(&GeometryInformation), "nrings"_a, "nph"_a, "ofs"_a, "stride"_a, "phi0"_a, "theta"_a, "wgt"_a )
    .def("nrings", &sharp_geom_info::nrings)
    .def("nph", &sharp_geom_info::nph, "iring"_a)
    .def("nphmax", &sharp_geom_info::nphmax)
    .def("theta", &sharp_geom_info::theta, "iring"_a)
    .def("cth", &sharp_geom_info::cth, "iring"_a)
    .def("sth", &sharp_geom_info::sth, "iring"_a)
    .def("phi0", &sharp_geom_info::phi0, "iring"_a)
    .def("pair", &sharp_geom_info::pair, "iring"_a)
    .def("clear_map", &sharp_geom_info::clear_map, "map"_a)
    .def("get_ring", &sharp_geom_info::get_ring, "weighted"_a, "iring"_a, "map"_a, "ringtmp"_a)
    .def("add_ring", &sharp_geom_info::add_ring, "weighted"_a, "iring"_a, "ringtmp"_a, "map"_a)
    .def("map2alm", &map2alm_ginfo, "map"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("alm2map", &alm2map_ginfo, "alm"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("map2alm_spin", &map2alm_spin_ginfo, "map"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
    .def("alm2map_spin", &alm2map_spin_ginfo, "alm"_a, "spin"_a, "lmax"_a, "mmax"_a, "nthreads"_a, "zbounds"_a)
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
  geom : Geometry
    A Scarf geometry following the HEALPix scheme
  )pbdoc", "nside"_a, "stride"_a);
}
