Scarf tutorial
===============

Prerequisites
--------------

This tutorial assumes that the reader is already familiar with other SHT libraries
such as `healpy <https://healpy.readthedocs.io/en/latest/>`_ or `ducc <https://gitlab.mpcdf.mpg.de/mtr/ducc/-/tree/ducc0>`_.


Basic map and alm conventions
-----------------------------

Scarf's alms and maps follow a similar convention to healpy's.

Alm arrays
~~~~~~~~~~

For spin-0 transforms, the alm array is assumed to be a one-dimensional array sorted by l, then by m (:math:`a_{0,0}, a_{1, 0} \ldots a_{l_{max}, 0}, a_{0, 1} \ldots a_{l_{max}, m_{max}}`)
As a result, users can use all `healpy.spht.Alm <https://healpy.readthedocs.io/en/latest/generated/healpy.sphtfunc.Alm.html>`_
functions on their scarf-created alms.

For higher spin transforms, the alm array is assumed to be a two-dimensional array of shape (2, hp.Alm.getsize(lmax, mmax)).
For spin-2 transforms, the first row would be the E alm array, and the second row would be the B alm array.

Map arrays
~~~~~~~~~~

For spin-0 transforms, the map array is assumed to be a one-dimensional array of length ``npix``, ``npix`` being the number of pixels
of the map.

For higher spin transforms, the map array is assumed to be a two-dimensional array of shape (2, npix).
For spin-2 functions, the first row would be the Q map, and the second row would be the U map.


HEALPix transforms
------------------

To perform transforms on HEALPix grid based maps, use the following functions:

- ``scarf.map2alm``
- ``scarf.map2alm_spin``
- ``scarf.alm2map``
- ``scarf.alm2map_spin``

Here is a minimal working example:

.. code-block:: python

   import scarf
   import numpy as np
   nside_mwe = 1
   npix_mwe = 12 * nside_mwe ** 2
   map_mwe = np.random.random(npix_mwe)
   lmax_mwe = 2
   
   scarf_alm = scarf.map2alm(
       map=map_mwe,
       nside=nside_mwe,
       lmax=lmax_mwe,
       mmax=lmax_mwe,
       nthreads=1,
       zbounds=[0, 1]
   )

The arguments of the functions are similar to healpy's with two notable additions:

- The ``nthreads`` parameter allow the user to specify the number of threads to use for the calculation (default = 1);
- The ``zbounds`` parameter allow the user to specify which portion of the map should be taken into account by the transform.
  Only rings whose latitude :math:`\theta` is in the interval ``zbounds`` = :math:`[\cos(\theta_1), \cos(\theta_2)]` are taken into account. 
   Setting ``zbounds = [0,1]`` thus restricts ``map2alm()`` to the northern hemisphere.


Transforms on custom geometries
-------------------------------

A scarf geometry is the specification of a map pixelization scheme.

Some examples are the HEALPix geometry,
the Gauss-Legendre (GL) geometry, or the equiangular geometry.

The user can create custom geometries and specify the number
of pixels, the weight and the location of each ring,
as well as the position of each ring's data in the map array used.

.. note::
   The pixels inside the ring are always uniformly spaced.


Declare a new geometry using the constructor: 

.. code-block:: python

   import scarf
   import numpy as np

   geom = scarf.Geometry(
      nrings=3,
      nph=[4, 4, 4],
      ofs=[0, 4, 8],
      stride=1,
      phi0=[np.pi/4, 0, np.pi/4],
      theta=[np.arccos(1-1/3), np.pi/2, np.arccos(1/3-1)],
      wgt=[np.pi/3]*3
   )

A geometry's attributes can be modified after the initial declaration. 

.. code-block:: python

   geom.theta[1] = [np.pi/3]

.. note::
   When creating a geometry with the constructor, the ring's placement are sorted by sin(``theta``) for performance reasons.
   This means the ordering of the input arrays may differ from the ordering of the attribute array.

Transforms are then available as methods:

.. code-block:: python

   from healpy import Alm

   mmax = 2
   lmax = 2

   scarf_alm = geom.alm2map(
       alm=np.random.random(Alm.getsize(lmax, mmax)),
       lmax=lmax,
       mmax=lmax
   )

Phase transforms
~~~~~~~~~~~~~~~~

A phase is the DUCC hidden representation of Legendre coefficients inside transforms.
If we define the spherical harmonics as

.. math::

   Y_{l m}(\theta, \phi) := P_{lm}(\cos \theta) e^{i m\phi}

with :math:`P_{lm}` being the associated Legendre polynomial, we can define a phase array:

.. math::

  \forall 1 \le j \le N_{rings}, 0 \le m \le m_{max},  \quad phase_j^m := \sum_{l=m}^{m_{max}} a_{lm} P_{lm}(\cos \theta_j)

Computing these values is a necessary step in both forward and backward transforms.

Methods to geometries exist to perform to compute
the phase values, and return them as an array.

.. code-block:: python

   from healpy import Alm

   mmax = 2
   lmax = 2

   scarf_alm = geom.alm2phase_spin(
       map=np.random.random(size=(2, Alm.getsize(lmax, mmax)))
       spin=2,
       lmax=lmax,
       mmax=lmax
   )


.. graphviz::

   digraph foo {
      label="Relationship between phase, map and alm";
      m [label="map"];
      p [label="phase"];
      a [label="alm"];
      { rank = same; m; a};
      { rank = same; p};
      m -> a [label="map2alm"];
      a -> m [label="alm2map"];
      a -> p [label="alm2phase"];
      m -> p [label="map2phase"];
      p -> a [label="phase2alm"];
      p -> m [label="phase2map"];
   }

.. note::
   For forward and backward transforms, DUCC's original phase arrays work differently.
   We've found the content of the phase between the two transforms differ by nph*weight
   for each ring. Scarf normalizes the output of the phase functions by dividing the phase obtained
   from map2phase and phase2alm by nph*weight.
