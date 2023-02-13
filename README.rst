import_serial
=============

Requirements: ccp4-python from CCP4 8.0 or later (it includes CCTBX)

Example of usage:

.. code ::

   $ ccp4-python -m import_serial --hklin 116720-721.lst-asdf-scale.hkl --half-dataset 116720-721.lst-asdf-scale.hkl1 116720-721.lst-asdf-scale.hkl2 --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90
   $ ccp4-python -m import_serial --hklin 116720-721.lst-asdf-scale.hkl --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 20 --dmin 1.65 --project protein --dataset 01

Test data are available in a separate repository: https://github.com/MartinMalyMM/import_serial_test_data

List of all options:

.. code ::

   $ ccp4-python -m import_serial --help
   
   usage: __main__.py [-h] --hklin hklin_file [--half-dataset hkl1 hkl2]
                      [--project PROJECT] [--crystal CRYST] [--dataset DATASET]
                      --spacegroup spacegroup --cell cell_a cell_b cell_c
                      cell_alpha cell_beta cell_gamma [--dmax D_MAX]
                      [--dmin D_MIN] [--nbins N_BINS]
   
   Covert hkl file(s) from CrystFEL to mtz and calculate statistics
   
   optional arguments:
     -h, --help            show this help message and exit
     --hklin hklin_file, --HKLIN hklin_file
                           Specify one merged hkl file from CrystFEL
     --half-dataset hkl1 hkl2
                           Two half-data-set merge files from CrystFEL (usually
                           .hkl1 and .hkl2)
     --project PROJECT     Project name
     --crystal CRYST       crystal
     --dataset DATASET     dataset
     --spacegroup spacegroup
                           Specify space group
     --cell cell_a cell_b cell_c cell_alpha cell_beta cell_gamma
                           Specify unit cell parameters divided by spaces, e.g.
                           60 50 40 90 90 90
     --dmax D_MAX, --lowres D_MAX
                           Low-resolution cutoff
     --dmin D_MIN, --highres D_MIN
                           High-resolution cutoff
     --nbins N_BINS, --nshells N_BINS
                           Number of resolution bins


Installation:

.. code ::

   $ ccp4-python -m pip install pairef --no-deps --upgrade --user

Test:

.. code ::

   $ cd test
   $ ccp4-python -m pytest -vv -s

Developed by Martin Maly, University of Southampton, `martin.maly@soton.ac.uk <mailto:martin.maly@soton.ac.uk>`_
