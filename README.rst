import_serial
=============

Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4.

Requirements: ccp4-python from CCP4 8.0 or later (it includes CCTBX)

Example of usage:

.. code ::

   $ ccp4-python -m import_serial --hklin 116720-721.lst-asdf-scale.hkl --half-dataset 116720-721.lst-asdf-scale.hkl1 116720-721.lst-asdf-scale.hkl2 --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90
   $ ccp4-python -m import_serial --hklin 116720-721.lst-asdf-scale.hkl --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 20 --dmin 1.65 --project protein --dataset 01

List of all options:

.. code ::

   $ ccp4-python -m import_serial --help
   
   usage: ccp4-python -m import_serial [-h] --hklin HKLIN [--spacegroup SPACEGROUP]
                   [--cell a b c alpha beta gamma] [--cellfile CELLFILE]
                   [--reference REFERENCE] [--half-dataset HKL1 HKL2]
                   [--project PROJECT] [--crystal CRYST] [--dataset DATASET]
                   [--dmax D_MAX] [--dmin D_MIN] [--nbins N_BINS]

   Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4

   optional arguments:
     -h, --help            show this help message and exit
     --hklin HKLIN, --HKLIN HKLIN
                           Specify merged mtz file from xia2.ssx or merged hkl
                           file from CrystFEL
     --spacegroup SPACEGROUP
                           Space group
     --cell a b c alpha beta gamma
                           Unit cell parameters divided by spaces, e.g. 60 50 40
                           90 90 90
     --cellfile CELLFILE   Cell file from CrystFEL
     --reference REFERENCE, --ref REFERENCE, --pdb REFERENCE, --cif REFERENCE, --mmcif REFERENCE
                           Reference file (PDB, mmCIF or MTZ) to provide
                           spacegroup and unit cell
     --half-dataset HKL1 HKL2
                           CrystFEL only: two half-data-set merge files (usually
                           .hkl1 and .hkl2)
     --project PROJECT     Project name
     --crystal CRYST       Crystal name
     --dataset DATASET     Dataset name
     --dmax D_MAX, --lowres D_MAX
                           Low-resolution cutoff
     --dmin D_MIN, --highres D_MIN
                           High-resolution cutoff
     --nbins N_BINS, --nshells N_BINS
                           Number of resolution bins

Installation
------------

The latest code from GitHub:

.. code ::

   $ ccp4-python -m pip install https://github.com/MartinMalyMM/import_serial/archive/master.zip --no-deps --upgrade --user

or - from PyPI:

.. code ::

   $ ccp4-python -m pip install import_serial --no-deps --upgrade --user

Tests
-----

.. code ::

   $ cd test
   $ ccp4-python -m pytest -vv -s

Test data are available in a separate repository: https://github.com/MartinMalyMM/import_serial_test_data

Developed by Martin Maly, University of Southampton, `martin.maly@soton.ac.uk <mailto:martin.maly@soton.ac.uk>`_
