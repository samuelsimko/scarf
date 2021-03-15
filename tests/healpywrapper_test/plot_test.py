import sys
import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import ducc0

from healpy import Alm


def get_nthreads():
    """ Returns the number of available threads on a posix/win based system """
    if sys.platform == "win32":
        return (int)(os.environ["NUMBER_OF_PROCESSORS"])
    return (int)(os.popen("grep -c cores /proc/cpuinfo").read())


NTHREADS = get_nthreads()

# Test data
NSIDE = 512
test_map = hp.read_map("./haslam408_dsds_Remazeilles2014.fits", dtype=np.float64)
LMAX = 2 ** 5


# Healpy
healpy_alms = hp.sphtfunc.map2alm(test_map, lmax=LMAX, mmax=LMAX, iter=1)

# DUCC
job = ducc0.sht.sharpjob_d()
job.set_nthreads(NTHREADS)
job.set_triangular_alm_info(LMAX, LMAX)
job.set_healpix_geometry(NSIDE)
ducc_alms = job.map2alm(test_map)


def imshow_alms(alms, LMAX, title="", show=True):
    """
    Plots the alms as pixels of an image.

            Parameters:
                    alms (np.ndarray): An array containing the alms
                                       ordered by the HEALPix standard scheme
                    LMAX (int): The maximum l value of the alms
                    title (string): The title of the figure
                    show (bool): If true, shows the figure
    """

    alms_2d = [
        [np.real(alms[Alm.getidx(LMAX, l - 1, m)]) for m in range(0, l)]
        + [0.0] * (LMAX - l + 1)
        for l in range(1, LMAX + 2)
    ]
    alms_2d_imag = [
        [np.imag(alms[Alm.getidx(LMAX, l - 1, m)]) for m in range(0, l)]
        + [0.0] * (LMAX - l + 1)
        for l in range(1, LMAX + 2)
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2)
    alm_re = ax1.imshow(alms_2d)
    alm_im = ax2.imshow(alms_2d_imag)

    # Hide the top part of the image with white rectangles
    for x in range(0, LMAX + 1):
        rect = patches.Rectangle([x - 0.5, x - 1.5], LMAX, 1, color="white")
        rect2 = patches.Rectangle([x - 0.5, x - 1.5], LMAX, 1, color="white")
        ax1.add_patch(rect)
        ax2.add_patch(rect2)

    ax1.set_xlabel("m")
    ax1.set_ylabel("l")
    ax2.set_xlabel("m")
    ax2.set_ylabel("l")
    fig.suptitle(title)
    ax1.title.set_text("Real values")
    ax2.title.set_text("Imaginary values")
    fig.colorbar(alm_re, ax=ax1)
    fig.colorbar(alm_im, ax=ax2)

    if show:
        plt.show()


if __name__ == "__main__":
    imshow_alms(healpy_alms, LMAX, "Healpy map2alm's produced alms", show=False)
    imshow_alms(ducc_alms, LMAX, "DUCC map2alm's produced alms", show=False)
    imshow_alms(
        ducc_alms - healpy_alms, LMAX, "Difference between the alms", show=False
    )

    # Relative difference plot
    x_values = range(len(healpy_alms))
    plt.figure(4)
    plt.plot(
        x_values,
        np.absolute(healpy_alms - ducc_alms)
        / np.maximum(np.absolute(healpy_alms), np.absolute(ducc_alms)),
        "r",
        label="DUCC vs healpy",
    )
    plt.xlabel("alms")
    plt.ylabel("Relative difference (" + r"$\frac{|x-y|}{\max(|x|, |y|)}$" + ")")
    plt.legend(loc="best")
    plt.title("Relative difference of healpy's and DUCC's produced alms")
    plt.show()
