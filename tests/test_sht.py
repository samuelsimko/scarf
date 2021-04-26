from __future__ import absolute_import, print_function

import healpy as hp
import scarf
import warnings

import numpy as np

warnings.filterwarnings(action="ignore", module="healpy")


def test_map2alm():
    nside = 128
    lmax = nside
    m = np.random.random(12 * nside ** 2)
    hp_alm = hp.sphtfunc.map2alm(m, lmax, iter=0, verbose=False)
    scarf_alm = scarf.map2alm(m, lmax, lmax, 1, [-1, 1])
    assert np.linalg.norm(hp_alm - scarf_alm) < 1e-7


def test_alm2map():
    nside = 128
    lmax = nside
    alm = np.random.random(hp.Alm.getsize(lmax)).astype(np.complex)
    hp_map = hp.sphtfunc.alm2map(alm, nside, lmax, verbose=False)
    scarf_map = scarf.alm2map(alm, nside, lmax, lmax, 1, [-1, 1])
    assert np.linalg.norm(hp_map - scarf_map) < 1e-7


def test_alm2map_spin():
    nside = 128
    lmax = nside
    almt = np.random.random(hp.Alm.getsize(lmax)).astype(np.complex)
    alme = np.random.random(hp.Alm.getsize(lmax)).astype(np.complex)
    almb = np.random.random(hp.Alm.getsize(lmax)).astype(np.complex)

    hp_map = hp.sphtfunc.alm2map(
        [almt, alme, almb], nside, lmax, pol=True, verbose=False
    )
    scarf_t = scarf.alm2map(almt, nside, lmax, lmax, 1, [-1, 1])
    [scarf_q, scarf_u] = scarf.alm2map_spin(
        [alme, almb], 2, nside, lmax, lmax, 1, [-1, 1]
    )
    assert np.linalg.norm(hp_map - [scarf_t, scarf_q, scarf_u]) < 1e-7


def test_map2alm_spin():
    nside = 128
    lmax = nside
    mmax = lmax
    mt = np.random.random(12 * nside ** 2)
    mq = np.random.random(12 * nside ** 2)
    mu = np.random.random(12 * nside ** 2)

<<<<<<< HEAD
    hp_map = hp.map2alm([mt, mq, mu], lmax, mmax, pol=True, verbose=False)
    scarf_map = scarf.map2alm_spin(np.array([mq, mq]), 0, nside, lmax, mmax, 1, [-1, 1])
    print(np.linalg.norm(hp_map - scarf_map))
    assert np.linalg.norm(hp_map - scarf_map) < 1e-7
=======
    hp_alm = hp.map2alm([mt, mq, mu], lmax, mmax, pol=True, verbose=False, iter=0)
    scarf_t = scarf.map2alm(mt, lmax, lmax, 1, [-1, 1])
    [scarf_e, scarf_b] = scarf.map2alm_spin(
        np.array([mq, mu]), 2, lmax, mmax, 1, [-1, 1]
    )
    assert np.linalg.norm(hp_alm - [scarf_t, scarf_e, scarf_b]) < 1e-7
>>>>>>> 815d2d0369f38c2fa4c593e589b2eb5eb2d6d589


def test_zbounds():
    nside = 128
    lmax = nside
    m = np.random.random(12 * nside ** 2)
    hp_alm = hp.sphtfunc.map2alm(m, nside, lmax, iter=0)
    hp_map = hp.sphtfunc.alm2map(hp_alm, nside, lmax)

    scarf_alm = scarf.map2alm(m, lmax, lmax, 1, [-1, 0])
    scarf_map = scarf.alm2map(scarf_alm, nside, lmax, lmax, 1, [-1, 0])
    assert (
        np.linalg.norm(hp_alm[6 * nside ** 2 : -1] - scarf_alm[6 * nside ** 2 : -1])
        < 1e-7
    )
    assert (
        np.linalg.norm(hp_alm[0 : 6 * nside ** 2] - scarf_alm[0 : 6 * nside ** 2]) > 1
    )


def test_geometry():
    nside = 1
    lmax = 4

    g = scarf.Geometry(
        nrings=3,
        nph=[4, 4, 4],
        ofs=[0, 4, 8],
        stride=1,
        phi0=[np.pi / 4, 0, np.pi / 4],
        theta=[np.arccos(1 - 1 / 3), np.pi / 2, np.arccos(1 / 3 - 1)],
        wgt=[np.pi / 3] * 3,
    )
    hg = scarf.healpix_geometry(1, 1)

    m = np.random.random(12 * nside ** 2)
    alm = np.random.random(hp.Alm.getsize(lmax)).astype(np.complex)

    g_alm = g.map2alm(m, lmax, lmax, 1, [-1, 1])
    hg_alm = hg.map2alm(m, lmax, lmax, 1, [-1, 1])

    g_map = g.alm2map(alm, lmax, lmax, 1, [-1, 1])
    hg_map = hg.alm2map(alm, lmax, lmax, 1, [-1, 1])
    assert np.linalg.norm(g_alm - hg_alm) < 1e-7
    assert np.linalg.norm(g_map - hg_map) < 1e-7


test_alm2map_spin()
