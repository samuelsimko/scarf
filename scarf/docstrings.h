/*
 * Docstrings for Scarf functions
 *
 *    Created 2021 by Samuel Simko.
 *    Based on the DUCC source code by Martin Reinecke.
 *
 * This file is subject to the terms and conditions of the MIT License.
 */

// #include <iostream>
// using namespace std;

const char * module_ds = "Spherical harmonics transform library for CMB lensing";

const char * map2alm_ds = R"(
  Computes alms from a temperature map.

  Parameters
  ----------
  map : np.array, shape (:math:`N_{pix}`,)
    The temperature map
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{alm}`)
    Temperature alm
  )";

const char * alm2map_ds = R"(
  Computes a Healpix temperature map from alm.

  Parameters
  ----------
  alm : np.array, shape (:math:`N_{alm}`,)
    The temperature alm.
  nside : int, scalar
    The number of sides for the output map. This is the number of pixels along the diagonal of a healpix-grid-base-pixel.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{pix}`)
    Temperature map
  )";

const char * map2alm_spin_ds = R"(
  Computes alms from a polarization map.

  Parameters
  ----------
  map : np.array, shape (:math:`N_{pix}`,)
    The polarisation map.
  spin : int, scalar
    The spin of the field.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{alm}`)
    Polarisation alm
    )";

const char * alm2map_spin_ds = R"(
  Computes a Healpix polarisation map from alm.

  Parameters
  ----------
  alm : np.array, shape (:math:`N_{alm}`,)
    The polarisation alm
  spin : int, scalar
    The spin of the field.
  nside : int, scalar
    The number of sides for the output map. This is the number of pixels along the diagonal of a healpix-grid-base-pixel.
  lmax: int, scalar
    The maximum angular momentum quantum number (multipole) of the powerspectrum.
  mmax: int, scalar
    The maximum magnetic quantum number of the powerspectrum.
  nthreads: int, scalar
    The number of threads for the computation.
  zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
    The latitudinal bounds. Only rings within zbounds will be calculated.

  Returns
  -------
  np.array, shape (:math:`N_{pix}`)
    Temperature map
    )";

const char * GL_wg_ds = R"(
  Computes Gauss-Legendre quadrature weights.

  Parameters
  ----------
  n: int, scalar
    The number of GL sample points used (integrates exactly degree :math:`2n - 1` polynomials; typically :math:`n=\ell_{\rm max} + 1`)

  Returns
  -------
  np.array, shape (:math:`n`)
  )";

const char * GL_xg_ds = R"(
  Computes Gauss-Legendre quadrature sample points. Output is ordered from 1 to -1.

  Parameters
  ----------
  n: int, scalar
    The number of GL sample points used (integrates exactly degree :math:`2n - 1` polynomials; typically :math:`n=\ell_{\rm max} + 1`)

  Returns
  -------
  np.array, shape (:math:`n`)
)";

const char * geometry_ds = R"(
  Creates a geometry specified by the user.

  Parameters
  ----------
  nrings : int, scalar
    The number of rings
  nph : int,  shape (nrings)
    TBD
  ofs : int,  shape (nrings)
    TBD
  stride : int, scalar
    The stride between two consecutive pixels on the map
  phi0 : float, shape (nrings)
    TBD
  theta : float,  shape (nrings)
    The latitude angle of each ring
  wgt : float,  shape (nrings)
    The weighting for each ring
    
  Attributes
  ----------
  theta : float,  shape (nrings)
    The latitude angle of each ring
  nph : int,  shape (nrings)
    TBD
  ofs : int,  shape (nrings)
    TBD
  stride : int, scalar
    The stride between two consecutive pixels on the map
  phi0 : float, shape (nrings)
    TBD
  wgt : float,  shape (nrings)
)";
const char * get_nrings_ds = R"(
    Returns the number of rings of the geometry.

    Parameters
    ----------

    Returns
    -------
    int
      The number of rings of the geometry
)";

const char * get_nph_ds = R"(
    Returns the number of pixels in the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    int
      The number of pixels in the ring.
      
)";
const char * get_nphmax_ds = R"(
    Returns the maximum number of pixels of a ring.

    Parameters
    ----------

    Returns
    -------
    int
      The maximum number of pixels of a ring.
)";
const char * get_theta_ds = R"(
    Returns the lattitude, in radiants, of a specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The latitude of the ring
)";
const char * get_phi0_ds = R"(
    Returns the longitude of the first pixel of a ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The longitude of the first pixel of the ring.
      
)";
const char * get_weight_ds = R"(
    Returns the weight of the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    double
      The weight of the ring. 
)";
const char * get_ofs_ds = R"(
    Returns the pixel offset of the specified ring.

    Parameters
    ----------
    iring: int
      The identifier of the ring.

    Returns
    -------
    int
      The pixel offset of the ring. 
)";
    
const char * map2alm_ginfo_ds = R"(
    Computes a temperature map given the alm.

    Parameters
    ----------
    alm : np.array, shape (:math:`N_{alm}`,)
      The temperature alm
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be taken into account

    Returns
    -------
    np.array, shape (:math:`N_{pix}`)
      Temperature map
)";
const char * alm2map_ginfo_ds = R"(
    Computes the alm of a temperature map.

    Parameters
    ----------
    map : np.array, shape (:math:`N_{pix}`,)
      The temperature map
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (:math:`N_{alm}`)
      Temperature map
)";
const char * map2alm_spin_ginfo_ds = R"(
    Computes the alm of a polarization map.

    Parameters
    ----------
    map : np.array, shape (2, :math:`N_{pix}`,)
      The polarization map
    spin : int, scalar
      The spin of the field.
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be taken into account.

    Returns
    -------
    np.array, shape (2, :math:`N_{alm}`)
      Polarization alm
)";
const char * alm2map_spin_ginfo_ds = R"(
    Computes a polarization map given the alm.

    Parameters
    ----------
    alm : np.array, shape (2, :math:`N_{alm}`,)
      The polarization alm
    spin : int, scalar
      The spin of the field.
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (2, :math:`N_{pix}`)
      Polarization map
)";


const char * alm2phase_ds = R"(
    Computes a temperature phase given the alm.

    Parameters
    ----------
    alm : np.array, shape (:math:`N_{alm}`)
      The temperature alm
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (1, 2*npairs, mmax+1)
      Temperature phase
)";

const char * phase2alm_ds = R"(
    Computes a temperature alm given the phase.

    Parameters
    ----------
    phase : np.array, shape (1, 2*npairs, mmax+1)
      The temperature phase
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (2, :math:`N_{alm}`)
      Temperature alm
)";

const char * phase2map_ds = R"(
    Computes a temperature map given the phase.

    Parameters
    ----------
    phase : np.array, shape (1, 2*npairs, mmax+1)
      The temperature phase
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (2, :math:`N_{pix}`)
      Temperature map
)";

const char * map2phase_ds = R"(
    Computes a temperature phase given the map.

    Parameters
    ----------
    map : np.array, shape (:math:`N_{pix}`) or (2, :math: `N_{pix}`)
      The temperature alm
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (1, 2*npairs, mmax+1)
      Temperature phase
)";

const char * alm2phase_spin_ds = R"(
    Computes a polarization phase given the alm.

    Parameters
    ----------
    alm : np.array, shape (2, :math:`N_{alm}`)
      The polarization alm
    spin : int, scalar
      The spin of the field.
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (1, 2*npairs, mmax+1)
      Polarization phase
)";

const char * phase2alm_spin_ds = R"(
    Computes a polarization alm given the phase.

    Parameters
    ----------
    phase : np.array, shape (1, 2*npairs, mmax+1)
      The polarization phase
    spin : int, scalar
      The spin of the field.
    lmax: int, scalar
      The maximum angular momentum quantum number (multipole) of the powerspectrum.
    mmax: int, scalar
      The maximum magnetic quantum number of the powerspectrum.
    nthreads: int, scalar
      The number of threads for the computation.
    zbounds: List, shape (2), zbounds :math:`\in \{[-1,1]\}^2`
      The latitudinal bounds. Only rings within zbounds will be calculated.

    Returns
    -------
    np.array, shape (2, :math:`N_{alm}`)
      Polarization alm
)";

const char * healpix_geometry_ds = R"(
  Creates a HEALPix geometry given the nside.

  Parameters
  ----------
  nside : int, scalar
    The nside of the HEALPix geometry
  stride: int, scalar
    The stride between two consecutive pixels on the map

  Returns
  -------
  Geometry
    A Scarf geometry following the HEALPix scheme
  )pbdoc", "nside"_a, "stride"_a);
)";

const char * gauss_geometry_ds = R"(
  Creates a gauss geometry given nrings and nphi.

  Parameters
  ----------
  nrings : int, scalar
    The number of rings of the geometry
  nphi: int, scalar
    The number of pixels in each ring

  Returns
  -------
  Geometry
    A Scarf geometry following the gauss scheme
)";
