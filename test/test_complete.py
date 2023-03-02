import pytest
import os
import shutil
# import filecmp
# import warnings
import urllib.request  # Python 3
from helper import run, tmp_environ
from import_serial import __version__


@pytest.mark.parametrize(
"project, files, spacegroup, cell, others",
[("116720-721.lst-asdf-scale-no-half-dataset",
  {"hklin": "116720-721.lst-asdf-scale.hkl"},
   "P21", "39.4 78.5 48.0 90 97.94 90", "--nbins 10 --dmin 1.65"),
 ("116720-721.lst-asdf-scale",
  {"hklin": "116720-721.lst-asdf-scale.hkl", "hklin1": "116720-721.lst-asdf-scale.hkl1", "hklin2": "116720-721.lst-asdf-scale.hkl2"},
   "P21", "39.4 78.5 48.0 90 97.94 90", "--nbins 10 --dmin 1.65")],
ids=["116720-721.lst-asdf-scale", "116720-721.lst-asdf-scale-no-half-dataset"])
def test_complete(project, files, spacegroup, cell, others, tmp_environ):
    print("\nPerforming tests...")
    urllib_wget = urllib.request.urlretrieve

    print("   Downloading test data from https://github.com/MartinMalyMM/import_serial_test_data ...")
    url_prefix = "https://raw.githubusercontent.com/MartinMalyMM/import_serial_test_data/main/"
    for f in files.values():
        filename, log = urllib_wget(url_prefix + f, f)
        print("   .")
        assert os.path.isfile(f)
    hklin = files["hklin"]
    arguments = "--HKLIN " + hklin + \
        " --spacegroup " + spacegroup + \
        " --cell " + cell + \
        " --project " + project + \
        " " + others
    if project == "116720-721.lst-asdf-scale":
        hklin1 = files["hklin1"]
        hklin2 = files["hklin2"]
        arguments += " --half-dataset " + hklin1 + " " + hklin2
    print("   Executing process...")
    cp = run(arguments)
    assert cp.returncode == 0
    # assert not cp.stderr
    if cp.stderr:
        print("STDERR: " + cp.stderr)

    expected_stdout = """Conversion of CrystFEL HKL file to MTZ and import into CCP4

"""
    if project == "116720-721.lst-asdf-scale":
        expected_stdout += """Command line arguments:
--HKLIN 116720-721.lst-asdf-scale.hkl --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --project 116720-721.lst-asdf-scale --nbins 10 --dmin 1.65 --half-dataset 116720-721.lst-asdf-scale.hkl1 116720-721.lst-asdf-scale.hkl2"""
    elif project == "116720-721.lst-asdf-scale-no-half-dataset":
        expected_stdout += """Command line arguments:
--HKLIN 116720-721.lst-asdf-scale.hkl --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --project 116720-721.lst-asdf-scale-no-half-dataset --nbins 10 --dmin 1.65"""

    expected_stdout += """

Input parameters:
  hklin 116720-721.lst-asdf-scale.hkl
  spacegroup P21
  cell [39.4, 78.5, 48.0, 90.0, 97.94, 90.0]
"""
    if project == "116720-721.lst-asdf-scale":
        expected_stdout += """  half_dataset ['116720-721.lst-asdf-scale.hkl1', '116720-721.lst-asdf-scale.hkl2']
"""
    expected_stdout += """  project """ + project + """
  d_min 1.65
  n_bins 10


DATA STATISTICS:
================

Overall values:

#observed: 2663073
#unique: 68546
completeness = 100.0 %
multiplicity = 38.9
<I> = 185.03
<I/sigma(I)> = 2.81
"""
    if project == "116720-721.lst-asdf-scale":
        expected_stdout += """CC1/2 = 0.835
CC* = 0.954
Rsplit = 0.313

Binned values:

   d_max   d_min     #obs   #uniq   mult.   %comp      <I>  <I/sI>   cc1/2     cc* r_split
   78.50    3.56   256683    6888   37.30  100.00    961.1     4.0   0.738   0.922   0.306
    3.56    2.82   267271    6856   39.00  100.00    403.0     3.9   0.810   0.946   0.291
    2.82    2.46   273621    6851   39.90  100.00    174.8     3.7   0.823   0.950   0.286
    2.46    2.24   279486    6890   40.60  100.00    119.0     3.6   0.818   0.949   0.289
    2.24    2.08   270929    6828   39.70  100.00     77.7     3.3   0.810   0.946   0.310
    2.08    1.96   260818    6852   38.10  100.00     47.5     2.9   0.787   0.939   0.354
    1.96    1.86   265410    6877   38.60  100.00     27.1     2.3   0.328   0.702   0.504
    1.86    1.78   269067    6831   39.40  100.00     16.4     1.8   0.609   0.870   0.592
    1.78    1.71   263852    6819   38.70  100.00     11.1     1.4   0.505   0.819   0.780
    1.71    1.65   255936    6854   37.30  100.00      7.6     1.1   0.276   0.658   1.141

MTZ file created: 116720-721.lst-asdf-scale_dataset.mtz
"""
    elif project == "116720-721.lst-asdf-scale-no-half-dataset":
        expected_stdout += """
Binned values:

   d_max   d_min     #obs   #uniq   mult.   %comp      <I>  <I/sI>
   78.50    3.56   256683    6888   37.30  100.00    961.1     4.0
    3.56    2.82   267271    6856   39.00  100.00    403.0     3.9
    2.82    2.46   273621    6851   39.90  100.00    174.8     3.7
    2.46    2.24   279486    6890   40.60  100.00    119.0     3.6
    2.24    2.08   270929    6828   39.70  100.00     77.7     3.3
    2.08    1.96   260818    6852   38.10  100.00     47.5     2.9
    1.96    1.86   265410    6877   38.60  100.00     27.1     2.3
    1.86    1.78   269067    6831   39.40  100.00     16.4     1.8
    1.78    1.71   263852    6819   38.70  100.00     11.1     1.4
    1.71    1.65   255936    6854   37.30  100.00      7.6     1.1

MTZ file created: 116720-721.lst-asdf-scale-no-half-dataset_dataset.mtz
"""
    stdout_all = cp.stdout.splitlines(True)
    assert "".join(stdout_all) == expected_stdout

    files_list = \
        [hklin,
         "program.xml"]
    if project == "116720-721.lst-asdf-scale":
        files_list = files_list + \
            [hklin1,
             hklin2,
             "116720-721.lst-asdf-scale-no-half-dataset_dataset.mtz",
             "116720-721.lst-asdf-scale-no-half-dataset_dataset.json",]
    for f in files_list:
        assert os.path.isfile(f)
        os.remove(f)  # Clean up
        pass
    return
