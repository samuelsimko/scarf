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


def map2alm_spin(map, spin, nside, lmax, mmax, nthreads, zbounds):
    offset = healpywrappercpp.sphtfunc.offset(nside, 1, zbounds)
    npix = healpywrappercpp.sphtfunc.get_npix(nside, 1, zbounds)
    return healpywrappercpp.sphtfunc.map2alm_spin(
        map[0, offset : offset + npix],
        map[1, offset : offset + npix],
        spin,
        nside,
        lmax,
        mmax,
        nthreads,
        zbounds,
    )


def alm2map_spin(alm, spin, nside, lmax, mmax, nthreads, zbounds):
    map1 = np.zeros(shape=(12 * nside ** 2))
    map2 = np.zeros(shape=(12 * nside ** 2))
    offset = healpywrappercpp.sphtfunc.offset(nside, 1, zbounds)
    npix = healpywrappercpp.sphtfunc.get_npix(nside, 1, zbounds)
    healpywrappercpp.sphtfunc.alm2map_spin(
        alm,
        spin,
        nside,
        lmax,
        mmax,
        nthreads,
        zbounds,
        map1[offset : offset + npix],
        map2[offset : offset + npix],
    )
    return (map1, map2)
