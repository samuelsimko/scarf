import sys
import time
import random
import healpywrapper

import numpy as np
import matplotlib.pyplot as plt

from healpy import Alm

NTHREADS = 12

NSIDE = 512
LMAX = NSIDE * 3
alm_array_size = Alm.getsize(LMAX)

if sys.argv[1] == "a2m":
    random_alm = np.random.random(alm_array_size) + 1j * np.random.random(
        alm_array_size
    )
elif sys.argv[1] == "m2a":
    random_map = np.random.random(12 * NSIDE ** 2)

NB_LOOP = 1
NB_POINTS = 200

steps = list(enumerate(np.linspace(-0.99, 1, NB_POINTS)))

times = [None] * NB_POINTS
for i in steps:
    starting_time = time.time()
    for j in range(NB_LOOP):
        if sys.argv[1] == "a2m":
            healpywrapper.alm2map(
                random_alm, NSIDE, LMAX, LMAX, NTHREADS, [-1, i[1]]
            )
        elif sys.argv[1] == "m2a":
            healpywrapper.map2alm(
                random_map, NSIDE, LMAX, LMAX, NTHREADS, [-1, i[1]]
            )
    end_time = time.time() - starting_time
    times[i[0]] = end_time

plt.plot(np.array(steps)[:, 1], times, "or")
plt.xlabel(r"$\cos(\theta)$")
plt.ylabel("Time in seconds")
plt.title(
    "Time taken by the "
    + sys.argv[1]
    + " algorithm for zbounds = "
    + r"$[-1, \cos(\theta)]$"
)
plt.show()
