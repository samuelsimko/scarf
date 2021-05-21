import healpy as hp
import numpy as np
import scarf as scarf


def test_map2phase2alm():
    nside = 128
    lmax = nside
    mmax = nside

    m = np.array(range(1, 12 * nside ** 2 + 1))
    geom = scarf.healpix_geometry(nside, 1)

    alm = hp.sphtfunc.map2alm(m, lmax=lmax, mmax=mmax, iter=0)
    phase_from_map = geom.map2phase(m, lmax, mmax, 1, [-1, 1])
    print(phase_from_map.shape)
    alm_from_map = geom.phase2alm(phase_from_map, lmax, mmax, 1, [-1, 1])

    assert np.linalg.norm(alm_from_map - alm) < 1e-7


def test_map2phase2alm_spin():
    nside = 128
    lmax = nside
    mmax = nside

    m = np.random.random(size=(3, 12 * nside ** 2))
    geom = scarf.healpix_geometry(nside, 1)
    alm = hp.sphtfunc.map2alm(m, lmax=lmax, mmax=mmax, iter=0, pol=True)

    phase_from_map = geom.map2phase(m[1:2], lmax, mmax, 1, [-1, 1])
    print(phase_from_map)
    alm_from_map = geom.phase2alm_spin(phase_from_map, 2, lmax, mmax, 1, [-1, 1])
    assert np.linalg.norm(alm_from_map[0] - alm[1]) < 1e-7


def test_alm2phase2map():
    nside = 128
    lmax = nside
    mmax = nside

    alm = np.random.random(hp.Alm.getsize(lmax, mmax)) * (1j + 2)
    map = hp.sphtfunc.alm2map(alm, nside, lmax=lmax, mmax=mmax)
    geom = scarf.healpix_geometry(nside, 1)
    phase_from_alm = geom.alm2phase(alm, lmax, mmax, 1, [-1, 1])
    print(phase_from_alm.shape)
    map_from_alm = geom.phase2map(phase_from_alm, lmax, mmax, 10, [-1, 1])

    assert np.linalg.norm(map_from_alm - map) < 1e-7


def test_alm2phase2map_spin():
    nside = 128
    lmax = nside
    mmax = nside

    alm = np.random.random(size=(3, hp.Alm.getsize(lmax, mmax))) * (1j + 2)
    map = hp.sphtfunc.alm2map(alm, nside, lmax=lmax, mmax=mmax, pol=True)
    geom = scarf.healpix_geometry(nside, 1)
    phase_from_alm = geom.alm2phase_spin(alm[1:], 2, lmax, mmax, 10, [-1, 1])

    print(phase_from_alm.shape)
    # phase_from_alm = np.reshape(phase_from_alm, (2, 512, 129))
    half = phase_from_alm[:, :, 0]
    print("half: ", half.shape)

    # phase_from_alm = np.array([phase_from_alm[:, :, 0], phase_from_alm[:, :, 1]])
    print(phase_from_alm.shape)
    map_from_alm = geom.phase2map(phase_from_alm, lmax, mmax, 10, [-1, 1])

    print(map_from_alm)
    print(map)

    print(np.linalg.norm(map_from_alm[0] - map[1]))
    assert np.linalg.norm(map_from_alm[0] - map[1]) < 1e-7


def test_alm2phase2alm():

    nside = 128
    lmax = nside // 1
    mmax = nside // 1

    map = np.random.random(12 * nside ** 2)
    # alm = np.random.random(hp.Alm.getsize(lmax, mmax))*(1j+2)
    alm = hp.sphtfunc.map2alm(map, lmax, iter=0)
    alm_back = hp.sphtfunc.map2alm(hp.sphtfunc.alm2map(alm, nside, lmax), lmax, iter=0)
    geom = scarf.healpix_geometry(nside, 1)
    phase_from_alm = geom.alm2phase(alm, lmax, mmax, 10, [-1, 1])
    alm_from_phase = geom.phase2alm(phase_from_alm, lmax, mmax, 10, [-1, 1])

    assert np.linalg.norm(alm_from_phase - alm_back) < 1e-5


def test_map2phase2map():
    nph = 128
    ntheta = 128
    lmax = 63
    mmax = 63

    geom = scarf.gauss_geometry(nph, ntheta)

    map = geom.alm2map(
        np.random.random(hp.Alm.getsize(lmax, mmax)) * (2 + 1j), lmax, mmax, 1
    )

    phase_from_map = geom.map2phase(map, lmax, mmax, 10, [-1, 1])
    map_from_phase = geom.phase2map(phase_from_map, lmax, mmax, 10, [-1, 1])

    assert np.linalg.norm(map - map_from_phase) < 1e-7
