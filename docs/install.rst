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