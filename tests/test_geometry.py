from __future__ import absolute_import, print_function

import healpy as hp
import scarf
import warnings

import numpy as np

warnings.filterwarnings(action="ignore", module="healpy")


def test_geometry_constructor():
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


def test_attributes():
    nside = 512
    g = scarf.healpix_geometry(512, 1)
    assert (
        np.linalg.norm(
            g.theta - [g.get_theta(iring) for iring in range(0, g.get_nrings())]
        )
        < 1e-7
    )
    assert (
        np.linalg.norm(
            g.weight - [g.get_weight(iring) for iring in range(0, g.get_nrings())]
        )
        < 1e-7
    )
    assert (
        np.linalg.norm(
            g.phi0 - [g.get_phi0(iring) for iring in range(0, g.get_nrings())]
        )
        < 1e-7
    )
    assert (
        np.linalg.norm(g.cth - [g.get_cth(iring) for iring in range(0, g.get_nrings())])
        < 1e-7
    )
    assert (
        np.linalg.norm(g.sth - [g.get_sth(iring) for iring in range(0, g.get_nrings())])
        < 1e-7
    )
    assert (
        np.linalg.norm(g.nph - [g.get_nph(iring) for iring in range(0, g.get_nrings())])
        < 1e-7
    )
