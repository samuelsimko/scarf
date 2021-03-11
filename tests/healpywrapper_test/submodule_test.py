import time
import sys
import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import ducc0
import healpywrapper_test

NSIDE = 512


def get_nthreads():
    """ Returns the number of available threads on a posix/win based system """
    if sys.platform == "win32":
        return (int)(os.environ["NUMBER_OF_PROCESSORS"])
    return (int)(os.popen("grep -c cores /proc/cpuinfo").read())


NTHREADS = get_nthreads()

LMAX = 100
# Test data
# test_map = hp.read_map("./haslam408_dsds_Remazeilles2014.fits", dtype=np.float64)

# """ Bigger data
test_map = hp.read_map("./COM_CMB_IQU-smica_2048_R3.00_hm1.fits", dtype=np.float64)
# NSIDE = 2048
# LMAX = NSIDE * 2
# """


NSIDE = 2048
test_map = np.random.random(size=(50_331_648))
LMAX = NSIDE * 2


# Healpy
starting_time = time.time()
healpy_alms = hp.sphtfunc.map2alm(test_map, lmax=LMAX, mmax=LMAX, iter=1)
end_time = time.time() - starting_time
print("healpy's map2alm time: ", end_time, " sec.")

# DUCC
starting_time = time.time()

job = ducc0.sht.sharpjob_d()
job.set_nthreads(NTHREADS)
job.set_triangular_alm_info(LMAX, LMAX)
job.set_healpix_geometry(NSIDE)
ducc_alms = job.map2alm(test_map)

end_time = time.time() - starting_time
print("DUCC's map2alm time: ", end_time, " sec.")


# healpywrapper_test
starting_time = time.time()
my_alms = healpywrapper_test.my_submodule.my_map2alm(
    test_map, NSIDE, LMAX, LMAX, NTHREADS
)
end_time = time.time() - starting_time
print("my_map2alm's time: ", end_time, " sec.")


# Tests
print("Norm L2 of (my_alms - healpy_alms): ", np.linalg.norm(my_alms - healpy_alms))
print(
    "Biggest difference between two alms (out of {} values): ".format(len(my_alms)),
    max(abs(my_alms - healpy_alms)),
)

print("Norm L2 of (my_alms - ducc_alms): ", np.linalg.norm(my_alms - ducc_alms))
print(
    "Biggest difference between two alms (out of {} values): ".format(len(my_alms)),
    max(abs(my_alms - ducc_alms)),
)

# Graph
"""
x_values = range(len(my_alms))
plt.plot(x_values, healpy_alms, 'r', label="healpy")
plt.plot(x_values, ducc_alms, 'g', label="ducc")
plt.plot(x_values, my_alms, 'b', label="my_map2alm")
plt.legend(loc='best')
plt.show()
"""
