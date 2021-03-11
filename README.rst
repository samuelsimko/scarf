==================
Healpy-wrapper
==================

About
-----

This repository is the home of my semester project for the course "Applications Informatiques" of the University of Geneva.

Eventually, this repository could become a wrapper to healpy, specifically designed for CMB applications.

Prerequisites
-------------

This repository uses pybind11, healpy and DUCC.
DUCC has several mandatory dependencies, make sure they are installed: 
https://gitlab.mpcdf.mpg.de/mtr/ducc


To install the repository prerequisites:

.. code-block:: console

    $ pip3 install --user pybind11 healpy pytest
    $ pip3 install --user git+https://gitlab.mpcdf.mpg.de/mtr/ducc.git

    $ git clone https://github.com/pybind/pybind11.git
    $ cd pybind11
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j8
    $ sudo make install
    

You also need to install openssl-devel or libssl-dev depending on your distribution with your
preferred package manager.

.. code-block:: console

    $ pip3 install --user pybind11 healpy pytest
    $ pip3 install --user git+https://gitlab.mpcdf.mpg.de/mtr/ducc.git

    $ git clone https://github.com/pybind/pybind11.git
    $ cd pybind11
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j8
    $ sudo make install

Usage
-----

To create the test binder .so library:

.. code-block:: console

    $ c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) pybind_test.cpp -o pybind_test$(python3-config --extension-suffix)

To test if the binder works:

.. code-block:: console

    $ cd healpy-wrapper/tests
    $ python binding_test.py

The last command should write "Hello world!" to STDOUT.

To create the healpywrapper_test module:

.. code-block:: console

    $ cd /tests/healpywrapper_test
    $ source compile.sh

The compilation takes one or two minutes.
This creates a shared object which can be imported in Python.
You can then compare DUCC, healpy, and healpywrapper_test's implementations of the function map2alm on random data:

.. code-block:: console

    $ python submodule_test.py
