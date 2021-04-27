.. currentmodule:: scarf

Geometry
=========

.. autosummary::
   :toctree: generated/

   Geometry


``scarf`` allows users to create their own map geometry, which allows for 
transforms on custom pixelations.

The class ``Geometry`` is used to store the pixelation information.

Call the constructor to create a new geometry:

.. code-block:: python


   import scarf

   geom = scarf.Geometry(
      nrings=3,
      nph=[4, 4, 4],
      ofs=[0, 4, 8],
      stride=1,
      phi0=[np.pi/4, 0, np.pi/4],
      theta=[np.arccos(1-1/3), np.pi/2, np.arccos(1/3-1)],
      wgt=[np.pi/3]*3
   )

Parameters

   - ``nrings`` is the number of iso-latitude rings.

   - ``nph`` is the number of pixels in each of these rings.

   - ``ofs`` is the index of the first pixel of each ring on a map which follows this geometry.

   - ``stride`` is the stride between pixels

   - ``phi0`` is the longitude angle, in radiants, of the first pixel of the ring.  It's value is in the interval [0, 2π].

   - ``theta`` is the lattitude angle, in radiants, of each ring. It's value is the interval [0, π], with 0 being the lattitude of the north pole.

   - ``wgt`` are the weights associated to each ring.




In this example, ``geom`` is equivalent to a HEALPix geometry of parameter ``nside`` = 1.

You can also create HEALPix geometries directly:

.. code-block:: python

   healpix_geom = scarf.healpix_geometry(nside=1, 1)

The ``Geometry`` class has methods to a variety of SHT transforms

.. code-block:: python

   alm = geom.map2alm(map, lmax=2, mmax=1, nthreads=12, zbounds = [-1, 0])

