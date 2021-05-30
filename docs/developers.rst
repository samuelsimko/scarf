Informations for Scarf developers
=================================

Good to see you want to improve ``scarf`` and it's functionalities!

General
**********

Scarf uses `pybind11 <https://pybind11.readthedocs.io/en/stable/>`_ to bind c++ code to python. 
It uses the `DUCC <https://gitlab.mpcdf.mpg.de/mtr/ducc/-/tree/ducc0>`_ library's c++ code for the SHT functions.

Scarf's c++ code can be seen as an interface between DUCC's c++ code and the scarf python package.

All of scarf's functions are compiled from c++ code onto your machine, creating a shared resource (.so)
As a result,

.. code-block:: console

   pip install --user --editable .

won't be sufficient to save your edit, as the .so file needs to be updated.
To do so,

   * Save any changes made to the files,
   * ``pip install --user .``

Another faster way to test your changes is to compile ``scarf.cc`` directly from your terminal:

.. code-block:: console

   cd scarf/
   c++ -std=c++17 -O3 -shared -march=native -Wall -fPIC $(python3 -m pybind11 --includes) scarf.cc -o scarf$(python3-config --extension-suffix) -fvisibility=hidden -I ../ducc/src

This will generate an .so object in your current directory, which can be imported to python (be careful to import the .so and not your scarf installation
by changing the relative import settings, or by executing the file in the same directory as the .so)

Library Organization
*********************

The c++ source code is located in the scarf folder.

In ``scarf.cc``:
   - Add the bindings to your functions

In ``docstrings.h``:
   - Add the docstrings to your functions (preferably in Napoleon's ``numpy`` format)

In ``functions.cc`` or in a new file:
   - Add your functions


Working in a Scarf branch
*************************

Create a branch, add, commit and push by executing,

.. code-block:: console

   git checkout -b <branchname>
   git add <whatever file changed>
   git commit -m '<your commit message>'
   git push -u origin <branchname>

Here, it is assumed that your remote repo is named ``origin``. For any subsequent pushes, it is sufficient to do,

.. code-block:: console

   git push


Working with an alternative DUCC branch
***************************************

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

It is slightly more tedious when working with submodules, espcially when switching between branches of the submodule.


The most simple solution is to execute the following, whenever switching between scarf or ducc branches,

.. code-block:: console

   git submodule sync
   git submodule update --remote

This guarantees that the correct branch and commit of the submodule is accessed, and that the local submodule files are the latest.

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


Formatting with Black
**********************

The repository uses Black for python file formatting.

When a python file is pushed to scarf's Github, the formatter `black <https://github.com/psf/black>`_ runs a check to see
if the python files are correctly formatted using Github Actions.
You can manually check your files by installing and running black:

.. code-block:: console

   pip install black
   python -m black --check <python file>

You can also automatically reformat it to fit black's format (don't forget to manually check the file afterwards)

.. code-block:: console

   python -m black <python file>
