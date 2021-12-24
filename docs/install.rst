Installation procedure
================================

``scarf`` can be installed with the following commands::

      git clone --recurse-submodules https://github.com/samuelsimko/scarf
      cd scarf
      git submodule update --remote
      pip install --user .

It is recommended to install ``scarf`` in a virtual environment.



Check your installation by importing the package and calling a function.
For instance, call the ``map2alm()`` function:

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


.. note::
   ``zbounds`` is the parameter controlling the latitude of the rings which are transformed.
   ``zbound = cos(latitude)``, where latitude goes from Pi to 0 radian.
   Setting ``zbounds = [0,1]`` thus restricts ``map2alm()`` to the northern hemisphere.

You can also run tests in the **tests** directory using ``pytest``.

.. code-block:: console

   pip install -U pytest
   python3 -m pytest tests
