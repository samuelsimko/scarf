# %%
import healpy as hp
import numpy as np
import scarf as scarfcpp

nside = 16
lmax = nside // 1
mmax = nside // 1

m = np.array(range(1, 12 * nside ** 2 + 1))
geom = scarfcpp.healpix_geometry(nside, 1)

alm = hp.sphtfunc.map2alm(m, lmax=lmax, mmax=mmax, iter=0)

phase_from_map = geom.map2phase(m, lmax, mmax, 10, [-1, 1])
print("phase_from_map:", phase_from_map)

phase_from_alm = geom.alm2phase(alm, lmax, mmax, 10, [-1, 1])
print("phase from alm: ", phase_from_alm)

alm_from_map = geom.phase2alm(phase_from_map, lmax, mmax, 10, [-1, 1])
print("alm from map: (with phase intermediate)", alm_from_map)
print(
    "phase2alm(alm2phase(map))) / map2alm(map)",
    alm_from_map / geom.map2alm(m, lmax, mmax, 10, [-1, 1]),
)

map_from_alm = geom.phase2map(phase_from_alm, lmax, mmax, 10, [-1, 1])

print("map from alm: (with phase intermediate)", map_from_alm)
print(
    "phase2map(map2phase(alm))) / alm2map(alm)",
    map_from_alm / geom.alm2map(alm, lmax, mmax, 10, [-1, 1]),
)

print("np.pi/(3*nside**2) : ", np.pi / (3 * nside ** 2))

# %%
