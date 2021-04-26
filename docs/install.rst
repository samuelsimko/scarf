Installation procedure
================================

We distinguish between,

      * a user, who will work with the main branch or a release of ``scarf``,
      * a developer, who uses his own scarf branches, as well as DUCC branches


User
-----


``scarf`` can be installed with the following commands::

      git clone --recursive-submodules https://github.com/samuelsimko/scarf
      cd scarf
      git submodule update --remote
      pip install --user .

It is recommended to install ``scarf`` in a virtual environment.



Developer
----------

Skip this section if you don't plan on working on scarf or DUCC

General
**********

The Python binder creates and installs a scarf package onto your machine by creating a shared resource (.so), thus,

.. code-block:: console

   pip install --user --editable .

won't be sufficient to have it execute your edit, as one needs to update the .so file.
To do so,

   * Save any changes you have made to the code,
   * ``pip install --user .``


Scarf branch
**************

Simply create your own branch, add, commit and push by executing,

.. code-block:: console

   git checkout -b <branchname>
   git add <whatever file changed>
   git commit -m '<your commit message>'
   git push -u origin <branchname>

Here, it is assumed that your remote repo is named ``origin``. For any subsequent pushes, it is sufficient to do,

.. code-block:: console

   git push



DUCC branch
************

We will discuss the following steps.
To make ``scarf`` work with an alternative ``DUCC`` branch, first,

      1. create a ``DUCC`` branch,
      2. update the git submodule information.

For any changes applied to ``DUCC`` branch,

      3. save changes, add, commit, and push to the alternative ``DUCC`` branch,
      4. let ``scarf`` point to the latest commit of the alternative ``DUCC`` branch.


**1. create a ``DUCC`` branch**

Access the DUCC repository, by locally changing the directory to **ducc/**. Assuming you're inside ``scarf``,

.. code-block:: console

   cd ducc/

   
Simply create your own branch, add, commit and push by executing,

.. code-block:: console

   git checkout -b <branchname>
   git add <whatever file changed>
   git commit -m '<your commit message>'
   git push -u origin <branchname>

Here, it is assumed that your remote repo is named ``origin``. For any subsequent pushes, it is sufficient to do,

   .. code-block:: console
   
      git push


**2. update the git submodule information**

Add the branch information of your submodule to ``.gitmodules`` by adding a **branch** parameter with the correct <branchname>

.. code-block:: rst

      [submodule "ducc"]
            path = ducc
            url = https://github.com/NextGenCMB/ducc/
            branch = <branchname>


@Sebastian TBD - not sure, do we need git submodule add upstream? testit


**3. save changes, add, commit, and push to the alternative DUCC branch**

.. code-block:: console

   -save <whatever file changed>-
   git add <whatever file changed>
   git commit -m '<your commit message>'
   git push -u origin <branchname>


**4. let ``scarf`` point to the latest commit of the alternative DUCC branch**

Inside the scarf root directory,

.. code-block:: console

   git add ducc/
   git commit -m '<your commit message>'
   git push


Switching between branches
***************************

Without submodules, switching between branches is simple,

.. code-block:: console

   git checkout <branchname>

locally replaces the current files with the files from the <branchname> branch and you are good to go.

It is slightly more tedious when working with submodules, espcially when you would like to switch the branch of the submodule.


The most simple solution is to execute the following, whenever switching between scarf or ducc branches,

.. code-block:: console

   git submodule sync
   git submodule update --remote

This guarantees that the correct branch and commit of the submodule is accessed, and that the local submodule files are the latest.




Minimal Working Example
========================

Import ``scarf`` (and ``numpy``, for the following example). It provides almost all spherical harmonic transforms
like healpy and follows a similar naming convention.

For instance, to calculate the alm from a given map, call the ``map2alm()`` function,

.. code-block:: python

   import scarf
   import numpy as np
   nside_mwe = 1
   npix_mwe = 12 * nside_mwe ** 2
   map_mwe = np.random.random(npix_mwe)
   lmax_mwe = 2
   
   scarf_alm = scarf.map2alm(
       map = map_mwe,
       nside = nside_mwe,
       lmax = lmax_mwe,
       mmax = lmax_mwe,
       nthreads = 1,
       zbounds = [0, 1])


``zbounds`` is the parameter controlling the latitude of the rings which are transformed.
``zbound = cos(latitude)``, where latitude goes from Pi to 0 radian.
Setting ``zbounds = [0,1]`` thus restricts ``map2alm()`` to the northern hemisphere.



Testing
================================

A basic pytest is currently executed upon each pull-request and push, for each branch.
To manually test the code with the existing test directory **tests**, install ``pytest``,

.. code-block:: console

   $ pip install -U pytest

and execute in the root directory of the repository,

.. code-block:: console

   $ pytest tests

or,

.. code-block:: console

   $ python3 -m pytest tests