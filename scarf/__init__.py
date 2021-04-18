import numpy as np
from scarfcpp import Geometry, healpix_geometry

__version__ = "0.1.0"


def map2alm(map, nside, lmax, mmax, nthreads, zbounds):
    offset = scarfcpp.offset(nside, 1, zbounds)
    npix = scarfcpp.get_npix(nside, 1, zbounds)
    return scarfcpp.map2alm(
        map[offset : offset + npix], nside, lmax, mmax, nthreads, zbounds
    )


def alm2map(alm, nside, lmax, mmax, nthreads, zbounds):
    map = np.zeros(shape=(12 * nside ** 2))
    offset = scarfcpp.offset(nside, 1, zbounds)
    npix = scarfcpp.get_npix(nside, 1, zbounds)
    scarfcpp.alm2map(
        alm, nside, lmax, mmax, nthreads, zbounds, map[offset : offset + npix]
    )
    return map


def map2alm_spin(map, spin, nside, lmax, mmax, nthreads, zbounds):
    offset = scarfcpp.offset(nside, 1, zbounds)
    npix = scarfcpp.get_npix(nside, 1, zbounds)
    return scarfcpp.map2alm_spin(
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
    offset = scarfcpp.offset(nside, 1, zbounds)
    npix = scarfcpp.get_npix(nside, 1, zbounds)
    scarfcpp.alm2map_spin(
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
