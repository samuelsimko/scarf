==================
Scarf
==================

About
-----

Scarf is a Spherical Harmonic Transform library designed for CMB lensing applications.

The repository uses `DUCC <https://gitlab.mpcdf.mpg.de/mtr/ducc>`_ for the calculations.

Installation
------------

To install the package:

.. code-block:: console

   $ git clone --recurse-submodules https://github.com/samuelsimko/scarf
   $ cd scarf
   $ pip install --user .


Features
--------
- Forward and backward SHTs with arbitrary spin and sky cut parameters
- Custom creation of map geometries with user-specified parameters


Minimal Working Example
-----------------------

Import `scarf` (and `numpy`, for the following example). It provides almost all spherical harmonic transforms
like healpy and follows a similar naming convention.

For instance, to calculate the alm from a given map, call the `map2alm()` function,

.. code-block:: python

   import scarf
   import numpy
   nside_mwe = 1
   map_mwe = np.random.random(12 * nside_mwe ** 2)
   lmax_mwe = 2
   
   scarf_alm = scarf.map2alm(
       map = map_mwe,
       nside = nside_mwe,
       lmax = lmax_mwe,
       mmax = lmax_mwe,
       nthreads = 1,
       zbounds = [-1, 1])


Testing
--------

A basic pytest is currently executed upon each pull-request and push, for each branch.
To manually test the code with the existing test, install pytest,

.. code-block:: console

   $ pip install -U pytest

and execute pytest in the root directory of the repo,

.. code-block:: console

   $ pytest

or,

.. code-block:: console

   $ python3 -m pytest tests