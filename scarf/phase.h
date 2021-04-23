/*
 *  This file is part of libsharp2.
 *
 *  libsharp2 is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp2 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp2; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* libsharp2 is being developed at the Max-Planck-Institut fuer Astrophysik */

/*! \file sharp_internal.h
 *  Internally used functionality for the spherical transform library.
 *
 *  Copyright (C) 2006-2020 Max-Planck-Society
 *  \author Martin Reinecke \author Dag Sverre Seljebotn
 */

#include <complex>
#include <vector>
#include "ducc0/sharp/sharp.h"
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/mav.h"


namespace ducc0 {

namespace detail_sharp {

using std::complex;

class phase_job : public sharp_job {

  public:
    mav<complex<double>, 3> *phase;

  public:
    phase_job (sharp_jobtype type_,
        size_t spin_, const vector<any> &alm_, const vector<any> &map, mav<complex<double>, 3> &phase,
        const sharp_geom_info &geom_info, const sharp_alm_info &alm_info, size_t flags_, int nthreads_)
      : sharp_job::sharp_job(type_, spin_, alm_, map, geom_info, alm_info, flags_, nthreads_) {
      this->phase = &phase;
  }
    
DUCC0_NOINLINE void execute()// override
  {
  size_t lmax = ainfo.lmax(),
         mmax = ainfo.mmax();

  auto norm_l = (type==SHARP_ALM2MAP_DERIV1) ?
    detail_sht::YlmBase::get_d1norm (lmax) :
    detail_sht::YlmBase::get_norm (lmax, spin);

/* clear output arrays if requested */
   init_output();

  size_t nchunks, chunksize;
  get_chunk_info(ginfo.npairs(), (spin==0) ? 128 : 64,
                 nchunks,chunksize);
//FIXME: needs to be changed to "nm"
//
  // auto phase = mav<dcmplx,3>::build_noncritical({2*chunksize,mmax+1,nmaps()});
  detail_sht::YlmBase ylmbase(lmax,mmax,spin);
  detail_sht::SHT_mode mode = (type==SHARP_MAP2ALM) ? detail_sht::MAP2ALM : 
                             ((type==SHARP_ALM2MAP) ? detail_sht::ALM2MAP : detail_sht::ALM2MAP_DERIV1);
/* chunk loop */
  for (size_t chunk=0; chunk<nchunks; ++chunk)
    {
    size_t llim=chunk*chunksize, ulim=min(llim+chunksize,ginfo.npairs());
    vector<detail_sht::ringdata> rdata(ulim-llim);
    for (size_t i=0; i<ulim-llim; ++i)
      {
      double cth = ginfo.cth(ginfo.pair(i+llim).r1);
      double sth = ginfo.sth(ginfo.pair(i+llim).r1);
      auto mlim = detail_sht::get_mlim(lmax, spin, sth, cth);
      size_t idx = 2*i;
      size_t midx = 2*i+1;
      if (ginfo.pair(i+llim).r2==~size_t(0)) midx=idx;
      rdata[i] = { mlim, idx, midx, cth, sth };
      }

/* map->phase where necessary */
    // map2phase(mmax, llim, ulim, phase);

    ducc0::execDynamic(ainfo.nm(), nthreads, 1, [&](ducc0::Scheduler &sched)
      {
      detail_sht::Ylmgen ylmgen(ylmbase);
      auto almtmp = mav<dcmplx,2>({lmax+2, nalm()});

      while (auto rng=sched.getNext()) for(auto mi=rng.lo; mi<rng.hi; ++mi)
        {
/* alm->alm_tmp where necessary */
        alm2almtmp(mi, almtmp, norm_l);
        ylmgen.prepare(ainfo.mval(mi));
        ducc0::detail_sht::inner_loop(mode, almtmp, *phase, rdata, ylmgen, mi);

/* alm_tmp->alm where necessary */
        almtmp2alm(mi, almtmp, norm_l);
        }
      }); /* end of parallel region */

/* phase->map where necessary */
    // phase2map (mmax, llim, ulim, phase);
    } /* end of chunk loop */
  }

};

using namespace std;

template<typename T> void phase_execute (sharp_jobtype type, size_t spin, const vector<any> &alm,
  mav<complex<T>, 3> &phase,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads)
  {
  size_t npix = 0;
  for (size_t i = 0; i < geom_info.nrings(); ++i){;
    npix += geom_info.nph(i);
  }
  std::vector<double> dummy_map(npix, 0);
  phase_job job(type, spin, alm, {&dummy_map[0]}, phase, geom_info, alm_info, flags, nthreads);
  job.execute();
  }

template<typename T> void phase_execute_phase2map (phase_job &job, mav<complex<T>, 3> &phase,
    sharp_geom_info &geom_info, int mmax, int spin)
  {
   size_t nchunks, chunksize;
   get_chunk_info(geom_info.npairs(), (spin==0) ? 128 : 64, 
       nchunks,chunksize);
   for (size_t chunk=0; chunk<nchunks; ++chunk){
     size_t llim=chunk*chunksize, ulim=min(llim+chunksize,geom_info.npairs());
     //job.init_output();
     job.phase2map(mmax, llim, ulim, phase);
    }
  }

template<typename T> void phase_execute_map2phase (phase_job &job, mav<complex<T>, 3> &phase,
    sharp_geom_info &geom_info, int mmax, int spin)
  {
    size_t nchunks, chunksize;
    get_chunk_info(geom_info.npairs(), (spin==0) ? 128 : 64, 
        nchunks,chunksize);

    for (size_t chunk=0; chunk<nchunks; ++chunk){
      size_t llim=chunk*chunksize, ulim=min(llim+chunksize,geom_info.npairs());
      job.map2phase(mmax, llim, ulim, phase);
    }
  }

template<typename T> void sharp_alm2phase(const std::complex<T> *alm, mav<complex<T>, 3> &phase,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Y, 0, {alm}, phase, geom_info, alm_info, flags, nthreads);
  }
template<typename T> void sharp_alm2phase_adjoint(std::complex<T> *alm, mav<complex<T>, 3> &phase,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Y, 0, {alm}, {phase}, geom_info, alm_info, flags, nthreads);
  }
template<typename T> void sharp_alm2phase_spin(size_t spin, const std::complex<T> *alm1, const std::complex<T> *alm2, mav<complex<T>, 3> &phase1, mav<complex<T>, 3> &phase2,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Yt, spin, {alm1, alm2}, {phase1, phase2}, geom_info, alm_info, flags, nthreads);
  }
template<typename T> void sharp_alm2phase_spin_adjoint(size_t spin, std::complex<T> *alm1, std::complex<T> *alm2, mav<complex<T>, 3> &phase1, mav<complex<T>, 3> &phase2,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Yt, spin, {alm1, alm2}, {phase1, phase2}, geom_info, alm_info, flags, nthreads);
  }
template<typename T> void sharp_phase2alm(std::complex<T> *alm, mav<complex<double>,3> &phase,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Yt, 0, {alm}, phase, geom_info, alm_info, flags, nthreads);
  }
template<typename T> void sharp_phase2alm_spin(size_t spin, std::complex<T> *alm1, std::complex<T> *alm2, mav<complex<T>, 3> &phase1, mav<complex<T>, 3> &phase2,
  const sharp_geom_info &geom_info, const sharp_alm_info &alm_info,
  size_t flags, int nthreads=1)
  {
  phase_execute(SHARP_Yt, spin, {alm1, alm2}, {phase1, phase2}, geom_info, alm_info, flags, nthreads);
  }

}
}
