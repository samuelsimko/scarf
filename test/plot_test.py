import sys
import os
import healpywrapper as hpwp

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from healpy import Alm

# Test data
NSIDE = 512
LMAX = NSIDE
test_map = hp.read_map("../tests/healpywrapper_test/haslam408_dsds_Remazeilles2014.fits", dtype=np.float64)
# test_map = np.random.random(12*NSIDE**2)
healpy_alms = hp.sphtfunc.map2alm(test_map, lmax=LMAX, mmax=LMAX, iter=1)


TEST_VALUE = float(sys.argv[2])


if sys.argv[1] == "m2a":
    # Healpy
    healpy_alms = hp.sphtfunc.map2alm(test_map, lmax=LMAX, mmax=LMAX, iter=1)

    # shift map by correct index.
    wrapper_map = test_map
    start_index = hp.pixelfunc.ang2pix(NSIDE, np.arccos(TEST_VALUE), 0)
    """
    print(start_index)
    end_index = hp.pixelfunc.ang2pix(NSIDE, np.arccos(-0.5), 0)
    wrapper_map = np.roll(wrapper_map, -1*start_index)
    wrapper_map = wrapper_map[0: end_index - start_index]
    print("wrapper_map", len(wrapper_map))
    print("test_map", len(test_map))
    """

    # Healpywrapper
    wrapper_alms = hpwp.map2alm(test_map, NSIDE, LMAX, LMAX, 12, [-0.5, TEST_VALUE])
    print("wrapper: ", wrapper_alms)
    print("healpy_alms: ", healpy_alms)
    wrapper_map = hp.sphtfunc.alm2map(wrapper_alms, nside=NSIDE)
    healpy_map = hp.sphtfunc.alm2map(healpy_alms, nside=NSIDE)

elif sys.argv[1] == "a2m":
    # alm2map test (compute alms and convert back to maps)
    random_alm = np.random.random(Alm.getsize(LMAX))+1j*np.random.random(Alm.getsize(LMAX))
    healpy_alms = hp.sphtfunc.map2alm(test_map, lmax=LMAX, mmax=LMAX, iter=1)
    wrapper_alms = healpy_alms
    wrapper_map = hpwp.alm2map(wrapper_alms, NSIDE, LMAX, LMAX, 12, [-0.5, TEST_VALUE])
    healpy_map = hp.sphtfunc.alm2map(healpy_alms, nside=NSIDE)


x_values = range(12*NSIDE**2)
x_values_2 = range(len(wrapper_map))
plt.plot(x_values, healpy_map, 'r', x_values_2, wrapper_map, 'b', markevery=10**4)
plt.show()

diff_map = healpy_map - wrapper_map
rel_diff_map = np.absolute(healpy_map - wrapper_map) / np.maximum(np.absolute(healpy_map), np.absolute(wrapper_map))


hp.mollview(healpy_map, title="Healpy")
hp.mollview(wrapper_map, title="Healpywrapper")
hp.mollview(diff_map, title="Diff")
hp.mollview(rel_diff_map, title="Relative diff")

# plt.show()

plt.figure(10)
plt.plot(
    x_values,
    np.absolute(healpy_map - wrapper_map)
    / np.maximum(np.absolute(healpy_map), np.absolute(wrapper_map)),
    "r",
    label="healpywrapper vs healpy",
    markevery=10**5
)
plt.xlabel("map pixels")
plt.ylabel("Relative difference (" + r"$\frac{|x-y|}{\max(|x|, |y|)}$" + ")")
plt.legend(loc="best")
plt.title("Relative difference of healpy's and healpywrapper's produced maps")

plt.figure(11)
plt.plot(x_values, healpy_map, 'r', x_values, wrapper_map, 'b', markevery=10**4)

plt.show()
