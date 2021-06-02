/*
 * Phase functions declarations
 *
 *    Created 2021 by Samuel Simko.
 *    Based on the DUCC source code by Martin Reinecke.
 *
 * This file is subject to the terms and conditions of the MIT License.
 */

#include "ducc0/infra/communication.cc"
#include "ducc0/infra/system.cc"
#include "ducc0/infra/threading.cc"
#include "ducc0/infra/types.cc"
#include "ducc0/math/geom_utils.cc"
#include "ducc0/math/pointing.cc"
#include "ducc0/math/pointing.h"
#include "ducc0/math/space_filling.cc"
#include "ducc0/sht/sharp.h"
#include <complex>
#include <vector>

#include "ducc0/bindings/pybind_utils.h"
#include "ducc0/healpix/healpix_base.cc"
#include "ducc0/healpix/healpix_base.h"
#include "ducc0/healpix/healpix_tables.cc"
#include "ducc0/healpix/healpix_tables.h"
#include <complex>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/string_utils.cc"
#include "ducc0/infra/string_utils.h"
#include "ducc0/math/constants.h"
#include "ducc0/sht/sharp.cc"
#include "ducc0/sht/sht.cc"
#include "ducc0/sht/sht.h"

namespace ducc0 {

namespace detail_sharp {

using std::complex;

// static size_t new_nchunks_max = 1;

/// Get the chunk information needed to compute the whole map in one chunk
static void get_singular_chunk_info(size_t ndata, size_t nmult, size_t &nchunks,
                                    size_t &chunksize, size_t spin = 1) {
  // size_t chunksize_min = 10;
  chunksize = ndata;
  nchunks = 1;
  /*
  chunksize = (ndata+new_nchunks_max-1)/new_nchunks_max;
  if (chunksize>=chunksize_min) // use max number of chunks
    chunksize = ((chunksize+nmult-1)/nmult)*nmult;
  else // need to adjust chunksize and nchunks
    {
    nchunks = (ndata+chunksize_min-1)/chunksize_min;
    chunksize = (ndata+nchunks-1)/nchunks;
    if (nchunks>1)
      chunksize = ((chunksize+nmult-1)/nmult)*nmult;
    }
  nchunks = (ndata+chunksize-1)/chunksize;
  */
}

class phase_job : public sharp_job {

public:
  /// the phase array
  mav<complex<double>, 3> *phase;

  phase_job(sharp_jobtype type_, size_t spin_, const vector<any> &alm_,
            const vector<any> &map, mav<complex<double>, 3> &phase,
            const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
            size_t flags_, int nthreads_)
      : sharp_job::sharp_job(type_, spin_, alm_, map, geom_info, alm_info,
                             flags_, nthreads_) {
    this->phase = &phase;
  }

  DUCC0_NOINLINE void execute() // override
  {
    size_t lmax = ainfo.lmax(), mmax = ainfo.mmax();

    auto norm_l = (type == SHARP_ALM2MAP_DERIV1)
                      ? detail_sht::YlmBase::get_d1norm(lmax)
                      : detail_sht::YlmBase::get_norm(lmax, spin);

    /* clear output arrays if requested */
    init_output();

    size_t nchunks, chunksize;
    get_singular_chunk_info(ginfo.npairs(), (spin == 0) ? 128 : 64, nchunks,
                            chunksize);
    // FIXME: needs to be changed to "nm"
    //
    detail_sht::YlmBase ylmbase(lmax, mmax, spin);
    detail_sht::SHT_mode mode =
        (type == SHARP_MAP2ALM)
            ? detail_sht::MAP2ALM
            : ((type == SHARP_ALM2MAP) ? detail_sht::ALM2MAP
                                       : detail_sht::ALM2MAP_DERIV1);
    /* chunk loop */
    for (size_t chunk = 0; chunk < nchunks; ++chunk) {
      size_t llim = chunk * chunksize,
             ulim = min(llim + chunksize, ginfo.npairs());
      vector<detail_sht::ringdata> rdata(ulim - llim);
      for (size_t i = 0; i < ulim - llim; ++i) {
        double cth = ginfo.cth(ginfo.pair(i + llim).r1);
        double sth = ginfo.sth(ginfo.pair(i + llim).r1);
        auto mlim = detail_sht::get_mlim(lmax, spin, sth, cth);
        size_t idx = 2 * i;
        size_t midx = 2 * i + 1;
        if (ginfo.pair(i + llim).r2 == ~size_t(0))
          midx = idx;
        rdata[i] = {mlim, idx, midx, cth, sth};
      }

      ducc0::execDynamic(ainfo.nm(), nthreads, 1, [&](ducc0::Scheduler &sched) {
        detail_sht::Ylmgen ylmgen(ylmbase);
        auto almtmp = mav<dcmplx, 2>({lmax + 2, nalm()});

        while (auto rng = sched.getNext())
          for (auto mi = rng.lo; mi < rng.hi; ++mi) {
            /* alm->alm_tmp where necessary */
            alm2almtmp(mi, almtmp, norm_l);
            ylmgen.prepare(ainfo.mval(mi));
            ducc0::detail_sht::inner_loop(mode, almtmp, *phase, rdata, ylmgen,
                                          mi);

            /* alm_tmp->alm where necessary */
            almtmp2alm(mi, almtmp, norm_l);
          }
      }); /* end of parallel region */
    } /* end of chunk loop */
  }
};

using namespace std;

template <typename T>
void phase_execute(sharp_jobtype type, size_t spin, const vector<any> &alm,
                   mav<complex<T>, 3> &phase, const sharp_geom_info &geom_info,
                   const sharp_alm_info &alm_info, size_t flags, int nthreads) {
  size_t npix = 0;
  for (size_t i = 0; i < geom_info.nrings(); ++i) {
    ;
    npix += geom_info.nph(i);
  }
  vector<double> dummy_map(npix, 0);
  if (spin == 0) {
    phase_job job(type, spin, alm, {&dummy_map[0]}, phase, geom_info, alm_info,
                  flags, nthreads);
    job.execute();
  } else {
    vector<double> dummy_map_spin(npix, 0);
    phase_job job(type, spin, alm, {&dummy_map[0], &dummy_map_spin[0]}, phase,
                  geom_info, alm_info, flags, nthreads);
    job.execute();
  }
}

template <typename T>
void phase_execute_phase2map(phase_job &job, mav<complex<T>, 3> &phase,
                             sharp_geom_info &geom_info, int mmax, int spin) {
  size_t nchunks, chunksize;
  get_singular_chunk_info(geom_info.npairs(), (spin == 0) ? 128 : 64, nchunks,
                          chunksize, spin);

  for (size_t chunk = 0; chunk < nchunks; ++chunk) {
    size_t llim = chunk * chunksize,
           ulim = min(llim + chunksize, geom_info.npairs());
    job.init_output();
    job.phase2map(mmax, llim, ulim, phase);
  }
}

template <typename T>
void phase_execute_map2phase(phase_job &job, mav<complex<T>, 3> &phase,
                             sharp_geom_info &geom_info, int mmax, int spin) {
  size_t nchunks, chunksize;
  get_singular_chunk_info(geom_info.npairs(), (spin == 0) ? 128 : 64, nchunks,
                          chunksize, spin);

  for (size_t chunk = 0; chunk < nchunks; ++chunk) {
    size_t llim = chunk * chunksize,
           ulim = min(llim + chunksize, geom_info.npairs());
    job.map2phase(mmax, llim, ulim, phase);
  }

  double factor;
  for (size_t xi = 0; xi < phase.shape(0); ++xi) {
    for (size_t yi = 0; yi < phase.shape(1); ++yi) {
      factor = geom_info.nph(yi) * geom_info.weight(yi);
      for (size_t zi = 0; zi < phase.shape(2); ++zi) {
        phase.v(xi, yi, zi) =
            complex<double>(phase(xi, yi, zi).real() / factor,
                            phase(xi, yi, zi).imag() / factor);
      }
    }
  }
}

template <typename T>
void sharp_alm2phase(const std::complex<T> *alm, mav<complex<T>, 3> &phase,
                     const sharp_geom_info &geom_info,
                     const sharp_alm_info &alm_info, size_t flags,
                     int nthreads = 1) {
  phase_execute(SHARP_Y, 0, {alm}, phase, geom_info, alm_info, flags, nthreads);
}
template <typename T>
void sharp_alm2phase_adjoint(std::complex<T> *alm, mav<complex<T>, 3> &phase,
                             const sharp_geom_info &geom_info,
                             const sharp_alm_info &alm_info, size_t flags,
                             int nthreads = 1) {
  phase_execute(SHARP_Y, 0, {alm}, phase, geom_info, alm_info, flags, nthreads);
}
template <typename T>
void sharp_alm2phase_spin(size_t spin, const std::complex<T> *alm1,
                          const std::complex<T> *alm2,
                          mav<complex<T>, 3> &phase,
                          const sharp_geom_info &geom_info,
                          const sharp_alm_info &alm_info, size_t flags,
                          int nthreads = 1) {
  phase_execute(SHARP_Y, spin, {alm1, alm2}, phase, geom_info, alm_info, flags,
                nthreads);
}
template <typename T>
void sharp_alm2phase_spin_adjoint(size_t spin, std::complex<T> *alm1,
                                  std::complex<T> *alm2,
                                  mav<complex<T>, 3> &phase,
                                  const sharp_geom_info &geom_info,
                                  const sharp_alm_info &alm_info, size_t flags,
                                  int nthreads = 1) {
  phase_execute(SHARP_Yt, spin, {alm1, alm2}, phase, geom_info, alm_info, flags,
                nthreads);
}
template <typename T>
void sharp_phase2alm(std::complex<T> *alm, mav<complex<double>, 3> &phase,
                     const sharp_geom_info &geom_info,
                     const sharp_alm_info &alm_info, size_t flags,
                     int nthreads = 1) {

  double factor;

  auto new_phase = mav<dcmplx, 3>::build_noncritical(
      {phase.shape(0), phase.shape(1), phase.shape(2)});
  for (size_t xi = 0; xi < phase.shape(0); ++xi) {
    for (size_t yi = 0; yi < phase.shape(1); ++yi) {
      factor = geom_info.nph(yi) * geom_info.weight(yi);
      for (size_t zi = 0; zi < phase.shape(2); ++zi) {
        new_phase.v(xi, yi, zi) =
            complex<double>(phase(xi, yi, zi).real() * factor,
                            phase(xi, yi, zi).imag() * factor);
      }
    }
  }

  phase_execute(SHARP_Yt, 0, {alm}, new_phase, geom_info, alm_info, flags,
                nthreads);
}
template <typename T>
void sharp_phase2alm_spin(size_t spin, std::complex<T> *alm1,
                          std::complex<T> *alm2, mav<complex<T>, 3> &phase,
                          const sharp_geom_info &geom_info,
                          const sharp_alm_info &alm_info, size_t flags,
                          int nthreads = 1) {

  auto new_phase = mav<dcmplx, 3>::build_noncritical(
      {phase.shape(0), phase.shape(1), phase.shape(2)});
  double factor;
  for (size_t xi = 0; xi < phase.shape(0); ++xi) {
    for (size_t yi = 0; yi < phase.shape(1); ++yi) {
      factor = geom_info.nph(yi) * geom_info.weight(yi);
      for (size_t zi = 0; zi < phase.shape(2); ++zi) {
        new_phase.v(xi, yi, zi) =
            complex<double>(phase(xi, yi, zi).real() * factor,
                            phase(xi, yi, zi).imag() * factor);
      }
    }
  }
  phase_execute(SHARP_Yt, spin, {alm1, alm2}, new_phase, geom_info, alm_info,
                flags, nthreads);
}

} // namespace detail_sharp
} // namespace ducc0
