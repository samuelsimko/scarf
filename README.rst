.. image:: scarflogo.jpg
   :target: scarflogo.jpg
   :alt: Scarf logo


==================
Scarf
==================

.. image:: https://github.com/samuelsimko/scarf/workflows/CI/badge.svg
   :target: https://github.com/samuelsimko/scarf/actions

.. image:: https://readthedocs.org/projects/scarfcmb/badge/?version=latest
   :target: https://readthedocs.org/projects/scarfcmb/?badge=latest

.. image:: https://img.shields.io/badge/licence-MIT-brightgreen
   :target: https://github.com/samuelsimko/scarf/blob/master/LICENSE

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black


About
-----

Scarf is a Spherical Harmonic Transform library designed for CMB lensing applications.

The repository uses `DUCC <https://github.com/mreineck/ducc>`_ for the calculations.

Using binders to the DUCC source code, it provides functions that allow users
to easily do transforms on portions of the sky or custom maps.

Features
--------

- Custom creation of map geometries with user-specified parameters
- Forward and backward SHTs with sky cut and arbitrary spin parameters
  on custom or predefined geometries
- Transforms from and to Legendre coefficients (phase)

Installation
------------

.. code-block:: console

   $ git clone --recurse-submodules https://github.com/samuelsimko/scarf
   $ pip install -r requirements.txt
   $ cd scarf
   $ pip install --user .


Documentation
-------------

Scarf's documentation can be found `here <https://scarfcmb.readthedocs.io/en/latest/>`_.
