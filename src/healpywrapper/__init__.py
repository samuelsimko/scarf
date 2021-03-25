import healpywrappercpp
import numpy as np


def map2alm(map, nside, lmax, mmax, nthreads, zbounds):
    offset = healpywrappercpp.sphtfunc.offset(nside, 1, zbounds)
    npix = healpywrappercpp.sphtfunc.get_npix(nside, 1, zbounds)
    return healpywrappercpp.sphtfunc.map2alm(
        map[offset : offset + npix], nside, lmax, mmax, nthreads, zbounds
    )


def alm2map(alm, nside, lmax, mmax, nthreads, zbounds):
    map = np.zeros(shape=(12 * nside ** 2))
    offset = healpywrappercpp.sphtfunc.offset(nside, 1, zbounds)
    npix = healpywrappercpp.sphtfunc.get_npix(nside, 1, zbounds)
    healpywrappercpp.sphtfunc.alm2map(
        alm, nside, lmax, mmax, nthreads, zbounds, map[offset : offset + npix]
    )
    return map
