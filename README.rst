=========
Equi7Grid
=========

.. image:: https://travis-ci.org/TUW-GEO/Equi7Grid.svg?branch=master
    :target: https://travis-ci.org/TUW-GEO/Equi7Grid

.. image:: https://coveralls.io/repos/github/TUW-GEO/Equi7Grid/badge.svg?branch=master
    :target: https://coveralls.io/github/TUW-GEO/Equi7Grid?branch=master

.. image:: https://badge.fury.io/py/Equi7Grid.svg
    :target: https://badge.fury.io/py/Equi7Grid

.. image:: https://readthedocs.org/projects/equi7grid/badge/?version=latest
    :target: https://equi7grid.readthedocs.io/en/latest/?badge=latest

A python class for working with Equi7Grid - how to convert to - how to use the tiling system - etc.

It's a python package that handles the geometric and geographic operations of a gridded and tiled projection system.
It was designed for data cubes ingesting satellite imagery and builds the basis for the Equi7Grid (see https://github.com/TUW-GEO/Equi7Grid).

A detailed documentation on the Equi7Grid definition is at:

~/docs/doc_files/

Overlays for visualisation in Google Earth can be found here:

~/docs/doc_files/google_earth_overlays/

Citation
========

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1048530.svg
   :target: https://doi.org/10.5281/zenodo.1048530

If you use the software in a publication then please cite it using the Zenodo DOI.
Be aware that this badge links to the latest package version.

Please select your specific version at https://doi.org/10.5281/zenodo.1048530 to get the DOI of that version.
You should normally always use the DOI for the specific version of your record in citations.
This is to ensure that other researchers can access the exact research artefact you used for reproducibility.

You can find additional information regarding DOI versioning at http://help.zenodo.org/#versioning

Installation
============

This package should be installable through pip:

.. code::

    pip install Equi7Grid

Installs for numba, scripy and gdal are required from conda or conda-forge.

Contribute
==========

We are happy if you want to contribute. Please raise an issue explaining what
is missing or if you find a bug. We will also gladly accept pull requests
against our master branch for new features or bug fixes.

Development setup
-----------------

For Development we recommend a ``conda`` environment.

Example installation script
---------------------------

The following script will install miniconda and setup the environment on a UNIX
like system. Miniconda will be installed into ``$HOME/miniconda``.

.. code::

   wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
   bash miniconda.sh -b -p $HOME/miniconda
   export PATH="$HOME/miniconda/bin:$PATH"
   conda create -n equi7grid_env python=2.7 numpy scipy pip gdal pyproj numba
   source activate equi7grid_env
   pip install pytileproj


This script adds ``$HOME/miniconda/bin`` temporarily to the ``PATH`` to do this
permanently add ``export PATH="$HOME/miniconda/bin:$PATH"`` to your ``.bashrc``
or ``.zshrc``

The last line in the example activates the ``equi7grid_env`` environment.

After that you should be able to run:

.. code::

    python setup.py test

to run the test suite.

Guidelines
----------

If you want to contribute please follow these steps:

- Fork the Equi7Grid repository to your account
- Clone the repository
- make a new feature branch from the Equi7Grid master branch
- Add your feature
- Please include tests for your contributions in one of the test directories.
  We use py.test so a simple function called test_my_feature is enough
- submit a pull request to our master branch

Note
====

This project has been set up using PyScaffold 2.5.11. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.