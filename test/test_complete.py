import pytest
import os
import shutil
# import filecmp
# import warnings
import urllib.request  # Python 3
from helper import run, tmp_environ
from import_serial import __version__


@pytest.mark.parametrize(
"project, files, others",
[("116720-721.lst-asdf-scale-no-half-dataset",
  {"hklin": "116720-721.lst-asdf-scale.hkl"},
   "--spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 10 --dmin 1.65 --wavelength 1.1"),
 ("116720-721.lst-asdf-scale",
  {"hklin": "116720-721.lst-asdf-scale.hkl", "hklin1": "116720-721.lst-asdf-scale.hkl1", "hklin2": "116720-721.lst-asdf-scale.hkl2"},
   "--spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 10 --dmin 1.65 --wavelength 1.1"),
# ("116720-721.lst-asdf-scale-auto-half-dataset",
#  {"hklin": "116720-721.lst-asdf-scale.hkl"},
#   "--spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 10 --dmin 1.65 --wavelength 1.1"),
 ("dials-xia2-ssx",
  {"hklin": "merged.mtz"},
   "--nbins 20")],
ids=["116720-721.lst-asdf-scale", "116720-721.lst-asdf-scale-no-half-dataset", "dials-xia2-ssx"])
def test_complete(project, files, others, tmp_environ):
    print("\nPerforming tests...")
    urllib_wget = urllib.request.urlretrieve
    files_not_now_allowed = ["116720-721.lst-asdf-scale.hkl1", "116720-721.lst-asdf-scale.hkl2"]
    for f in files_not_now_allowed:
        if os.path.isfile(f):
            os.remove(f)

    print("   Downloading test data from https://github.com/MartinMalyMM/import_serial_test_data ...")
    url_prefix = "https://raw.githubusercontent.com/MartinMalyMM/import_serial_test_data/main/"
    for f in files.values():
        filename, log = urllib_wget(url_prefix + f, f)
        print("   .")
        assert os.path.isfile(f)
    hklin = files["hklin"]
    arguments = "--HKLIN " + hklin + \
        " --project " + project + \
        " " + others
    # " --spacegroup " + spacegroup + \
    # " --cell " + cell + \
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

    expected_stdout = """
Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4

"""
    if project == "116720-721.lst-asdf-scale":
        expected_stdout += """Command line arguments:
--HKLIN 116720-721.lst-asdf-scale.hkl --project 116720-721.lst-asdf-scale --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 10 --dmin 1.65 --wavelength 1.1 --half-dataset 116720-721.lst-asdf-scale.hkl1 116720-721.lst-asdf-scale.hkl2"""
    elif project == "116720-721.lst-asdf-scale-no-half-dataset":
        expected_stdout += """Command line arguments:
--HKLIN 116720-721.lst-asdf-scale.hkl --project 116720-721.lst-asdf-scale-no-half-dataset --spacegroup P21 --cell 39.4 78.5 48.0 90 97.94 90 --nbins 10 --dmin 1.65 --wavelength 1.1"""
    if "116720-721.lst-asdf-scale" in project:
        expected_stdout += f"""

Input parameters:
  hklin 116720-721.lst-asdf-scale.hkl
"""
        if project == "116720-721.lst-asdf-scale":
            expected_stdout += """  half_dataset ['116720-721.lst-asdf-scale.hkl1', '116720-721.lst-asdf-scale.hkl2']
"""
        expected_stdout += f"""  wavelength 1.1
  spacegroup P21
  cell [39.4, 78.5, 48.0, 90.0, 97.94, 90.0]
  d_min 1.65
  n_bins 10
  project {project}


DATA STATISTICS:
================

Overall values:

#observed: 2663073
#unique: 68546
completeness = 100.00 %
multiplicity = 38.85
<I> = 185.0
<I/sigma(I)> = 2.8
"""
        if project == "116720-721.lst-asdf-scale":
            expected_stdout += """CC1/2 = 0.835
CC* = 0.954
Rsplit = 0.313

Binned values:

   d_max   d_min     #obs   #uniq   mult.   %comp      <I>  <I/sI>   cc1/2     cc* r_split
   78.50    3.56   256683    6888   37.27  100.00    961.1     4.0   0.738   0.922   0.306
    3.56    2.82   267271    6856   38.98  100.00    403.0     3.9   0.810   0.946   0.291
    2.82    2.46   273621    6851   39.94  100.00    174.8     3.7   0.823   0.950   0.286
    2.46    2.24   279486    6890   40.56  100.00    119.0     3.6   0.818   0.949   0.289
    2.24    2.08   270929    6828   39.68  100.00     77.7     3.3   0.810   0.946   0.310
    2.08    1.96   260818    6852   38.06  100.00     47.5     2.9   0.787   0.939   0.354
    1.96    1.86   265410    6877   38.59  100.00     27.1     2.3   0.328   0.702   0.504
    1.86    1.78   269067    6831   39.39  100.00     16.4     1.8   0.609   0.870   0.592
    1.78    1.71   263852    6819   38.69  100.00     11.1     1.4   0.505   0.819   0.780
    1.71    1.65   255936    6854   37.34  100.00      7.6     1.1   0.276   0.658   1.141

MTZ file created: 116720-721.lst-asdf-scale_dataset.mtz
"""
        elif project == "116720-721.lst-asdf-scale-no-half-dataset":
            expected_stdout += """
Binned values:

   d_max   d_min     #obs   #uniq   mult.   %comp      <I>  <I/sI>
   78.50    3.56   256683    6888   37.27  100.00    961.1     4.0
    3.56    2.82   267271    6856   38.98  100.00    403.0     3.9
    2.82    2.46   273621    6851   39.94  100.00    174.8     3.7
    2.46    2.24   279486    6890   40.56  100.00    119.0     3.6
    2.24    2.08   270929    6828   39.68  100.00     77.7     3.3
    2.08    1.96   260818    6852   38.06  100.00     47.5     2.9
    1.96    1.86   265410    6877   38.59  100.00     27.1     2.3
    1.86    1.78   269067    6831   39.39  100.00     16.4     1.8
    1.78    1.71   263852    6819   38.69  100.00     11.1     1.4
    1.71    1.65   255936    6854   37.34  100.00      7.6     1.1

MTZ file created: 116720-721.lst-asdf-scale-no-half-dataset_dataset.mtz
"""
    elif project == "dials-xia2-ssx":
        expected_stdout += """Command line arguments:
--HKLIN merged.mtz --project dials-xia2-ssx --nbins 20

Input parameters:
  hklin merged.mtz
  n_bins 20
  project dials-xia2-ssx


DATA STATISTICS:
================

Overall values:

#observed: 342219
#unique: 35588
completeness = 87.30 %
multiplicity = 9.62
<I> = 107.2
<I/sigma(I)> = 11.5
CC1/2 = 0.945
CC* = 0.986
Rsplit = 0.281

Binned values:

   d_max   d_min     #obs   #uniq   mult.   %comp      <I>  <I/sI>   cc1/2     cc* r_split
"""
        expected_stdout_list = []
        expected_stdout_list.append("   68.16    4.30    40753    2169   18.79")
        expected_stdout_list.append("527.2    45.9   0.937   0.984   0.179")
        expected_stdout_list.append("    4.29    3.41    31332    2079   15.07")
        expected_stdout_list.append("484.1    42.3   0.933   0.982   0.177")
        expected_stdout_list.append("    3.41    2.98    28746    2055   13.99")
        expected_stdout_list.append("243.4    25.7   0.907   0.975   0.213")
        expected_stdout_list.append("    2.98    2.71    26812    2059   13.02")
        expected_stdout_list.append("134.6    16.7   0.869   0.964   0.248")
        expected_stdout_list.append("    2.71    2.52    25981    2032   12.79")
        expected_stdout_list.append("90.1    12.6   0.839   0.955   0.286")
        expected_stdout_list.append("    2.52    2.37    25295    2055   12.31")
        expected_stdout_list.append("68.6    10.3   0.834   0.954   0.301")
        expected_stdout_list.append("    2.37    2.25    24690    2041   12.10")
        expected_stdout_list.append("59.7     9.2   0.809   0.946   0.324")
        expected_stdout_list.append("    2.25    2.15    24789    2040   12.15")
        expected_stdout_list.append("49.0     8.3   0.805   0.944   0.351")
        expected_stdout_list.append("    2.15    2.07    23287    2004   11.62")
        expected_stdout_list.append("41.9     6.9   0.766   0.931   0.389")
        expected_stdout_list.append("    2.07    2.00    23238    2045   11.36")
        expected_stdout_list.append("32.4     5.9   0.687   0.903   0.463")
        expected_stdout_list.append("    2.00    1.93    18667    2022    9.23")
        expected_stdout_list.append("25.6     4.5   0.596   0.864   0.586")
        expected_stdout_list.append("    1.93    1.88    14318    1995    7.18")
        expected_stdout_list.append("21.2     3.4   0.376   0.739   0.819")
        expected_stdout_list.append("    1.88    1.83    10341    2016    5.13")
        expected_stdout_list.append("15.0     2.1   0.184   0.557   1.414")
        expected_stdout_list.append("    1.83    1.78     8009    1952    4.10")
        expected_stdout_list.append("11.5     1.5   0.088   0.403   1.892")
        expected_stdout_list.append("    1.78    1.74     6028    1921    3.14")
        expected_stdout_list.append("8.8     0.9   0.054   0.320   2.815")
        expected_stdout_list.append("    1.74    1.71     4402    1767    2.49")
        expected_stdout_list.append("7.3     0.7   0.044   0.290   4.701")
        expected_stdout_list.append("    1.71    1.67     2858    1470    1.94")
        expected_stdout_list.append("5.8     0.6   0.020   0.197   5.766")
        expected_stdout_list.append("    1.67    1.64     1779    1122    1.59")
        expected_stdout_list.append("5.2     0.5   0.069   0.360   4.897")
        expected_stdout_list.append("    1.64    1.61      725     586    1.24")
        expected_stdout_list.append("2.0     0.2   0.037   0.267 -37.419")
        expected_stdout_list.append("    1.61    1.58      169     158    1.07")
        expected_stdout_list.append("6.4     0.2  -0.523   0.000  -6.458")

    stdout_all = cp.stdout.splitlines(True)
    # print("".join(stdout_all))  # for debugging
    expected_stdout_list = expected_stdout.splitlines(keepends=True)
    for i, line in enumerate(expected_stdout_list):
        # print(i)            # for debugging
        assert expected_stdout_list[i] == stdout_all[i]
    # assert expected_stdout in "".join(stdout_all)
    if project == "dials-xia2-ssx":
        for halfline in expected_stdout_list:
            assert halfline in "".join(stdout_all)

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
