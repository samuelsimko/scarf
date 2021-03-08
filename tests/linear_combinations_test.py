#!/usr/bin/env python

import numpy as np
import healpy as hp
from healpy import Alm
import matplotlib.pyplot as plt
import ducc0

NSIDE = 2 ** 5
LMAX = 10

# Get spherical coordinates of RING convention
[theta, phi] = hp.pixelfunc.pix2ang(NSIDE, range(hp.pixelfunc.nside2npix(NSIDE)))

# Some of the first spherical harmonics
y_0_0 = lambda y: 1 / 2 * np.sqrt(1 / np.pi)
y_1_1 = lambda y: -1 / 2 * np.sqrt(3 / (2 * np.pi)) * np.sin(y[0]) * np.exp(1j * y[1])
y_2_m2 = (
    lambda y: 1
    / 4
    * np.sqrt(15 / (2 * np.pi))
    * (np.exp(-2j * y[1]) * np.sin(y[0]) ** 2)
)
y_2_0 = lambda y: 1 / 4 * np.sqrt(5 / np.pi) * (3 * np.cos(y[0]) ** 2 - 1)

# Custom function to mess around with
f = (
    lambda x: np.pi * y_1_1(x)
    + np.exp(1) * y_2_0(x)
    + (42 + 73j) * y_0_0(x)
    - 69 * y_2_m2(x)
)

# Transform f to HEALPix map
f_map = np.array([f([theta[i], phi[i]]) for i in range(len(theta))])

"""
# Visualise map using mollweide projection
hp.mollview(f_map)
plt.show()

# Compute and show power spectrum (code from the readthedocs)
cl = hp.sphtfunc.anafast(f_map, lmax=LMAX)
print(cl)

ell = np.arange(len(cl))
plt.figure(figsize=(10, 5))
plt.plot(ell, ell * (ell + 1) * cl)
plt.xlabel("$\ell$")
plt.ylabel("$\ell(\ell+1)C_{\ell}$")
plt.grid()
plt.show()
exit()
# plt.show()
"""

# Compute alms with healpy and DUCC
f_map_alm_healpy = hp.sphtfunc.map2alm(f_map, lmax=LMAX)

job = ducc0.sht.sharpjob_d()
job.set_triangular_alm_info(LMAX, LMAX)
nlon = 4096
nlat = LMAX + 1
job.set_healpix_geometry(NSIDE)
f_map_alm_ducc = job.map2alm(f_map)

"""
x_axis = range(len(f_map_alm_ducc))
plt.plot(x_axis, f_map_alm_ducc, 'y', x_axis, f_map_alm_healpy, 'r')
print(np.linalg.norm(f_map_alm_ducc - f_map_alm_healpy))
plt.show()
# The alms seem to be the same.
"""


f_map_healpy = hp.sphtfunc.alm2map(f_map_alm_healpy, nside=NSIDE)
f_map_ducc = hp.sphtfunc.alm2map(f_map_alm_ducc, nside=NSIDE)

x_axis = range(len(f_map_ducc))
fig, axs = plt.subplots(3)
fig.suptitle("Map after forward and backward transformation")
axs[2].plot(x_axis, f_map_ducc, "y", label="ducc")
axs[1].plot(x_axis, f_map_healpy, "r", label="healpy")
axs[0].plot(x_axis, f_map, "b", label="Original map")
axs[0].legend()
axs[1].legend()
axs[2].legend()
plt.show()

# Compute and show errors
print(np.linalg.norm(f_map_healpy - f_map))
print(np.linalg.norm(f_map_ducc - f_map))

# Did the coefficients stay the same ?
print("a_0_0: ", f_map_alm_healpy[Alm.getidx(LMAX, 0, 0)])
print("a_1_1: ", f_map_alm_healpy[Alm.getidx(LMAX, 1, 1)])
# (using conjugate formula to calculate a_2_m2)
print("a_2_m2: ", (-1) ** 2 * np.conj(f_map_alm_healpy[Alm.getidx(LMAX, 2, 2)]))
print("a_2_0: ", (f_map_alm_healpy[Alm.getidx(LMAX, 2, 0)]))
