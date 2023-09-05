# coding: utf-8
import argparse
import os
import sys
from pathlib import Path
import subprocess
import traceback
import math
import pandas as pd
from math import sqrt
import json
try:
    from cctbx import miller, crystal, uctbx, sgtbx, xray
    from cctbx.array_family import flex
    from iotbx import mtz, reflection_file_reader
except ModuleNotFoundError:
    print("WARNING: ModuleNotFoundError: Module CCTBX was not found.")
except ImportError:
    print("WARNING: ImportError: Module CCTBX was not found.")


def hkl_strip(hklin):
    """Keeps only lines in the format:
    int int int float whatever
    in the `hklin` file and saves them in the file `hklin_strip`.
    This routine checks only the lines in the beginning and end of
    the file."""
    with open(hklin, "r") as hklfile:
        lines = hklfile.readlines()
    line_start = 0
    line_end = len(lines)
    for i in range(len(lines)):
        try:
            int(lines[i].split()[0])
            int(lines[i].split()[1])
            int(lines[i].split()[2])
            float(lines[i].split()[3])
            break
        #except ValueError or IndexError:
        #except IndexError:
        #    sys.stderr.write(f"ERROR: Provided file {hklin} is in an incorrect format.\n"
        #                     "Aborting.\n")
        #    sys.exit(1)
        except:
            line_start = i + 1
    for i in reversed(range(len(lines))):
        try:
            int(lines[i].split()[0])
            int(lines[i].split()[1])
            int(lines[i].split()[2])
            float(lines[i].split()[3])
            break
        except:
            line_end = i
    lines_strip = lines[line_start:line_end]
    hkltmp = hklin + "_tmp"
    with open(hkltmp, "w") as hklfile:
        hklfile.writelines(lines_strip)
    return hkltmp


def which(program):
    """Checks if `program` exists and finds its location. Analogy of the
    `which` GNU/Linux command.
    Args:
        program (str): Name of an executable
    Returns:
        str: Path of an executable location
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program) or is_exe(program + ".exe"):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
            exe_file = os.path.join(path, program + ".exe")
            if is_exe(exe_file):
                return exe_file
    return None


def calc_rsplit(I1, I2):
    den = flex.sum(I1 + I2)
    if den == 0:
        return 0
    return math.sqrt(2) * flex.sum(flex.abs(I1 - I2)) / den


def calc_CCstar(CChalf):
    try:
        CCstar = sqrt(2 * float(CChalf) / (1 + float(CChalf)))
    except ValueError:
        CCstar = 0  # DIRTY HACK
    return CCstar


def get_miller_array_crystfel(hklin, cs, values="I", d_max=0, d_min=0):
    assert values == "I" or values == "nmeas"
    # read data from the text file
    # expected format of a fixed-width .hkl file:
    #    h    k    l          I    phase   sigma(I)   nmeas
    hklin_strip = hkl_strip(hklin)
    hklin_df = pd.read_csv(
        hklin_strip, header=None, index_col=False, sep=r'\s+',
        names=("h", "k", "l", "I", "phase", "sigma(I)", "nmeas"))
    h = flex.int(hklin_df["h"])
    k = flex.int(hklin_df["k"])
    l = flex.int(hklin_df["l"])
    indices = flex.miller_index(h,k,l)

    if values == "I":
        I = flex.double(hklin_df["I"])
        sig = flex.double(hklin_df["sigma(I)"])
        m = miller.array(
            miller_set=miller.set(cs, indices),
            data=I,
            sigmas=sig,
        )
        m.set_observation_type_xray_intensity()
    elif values == "nmeas":
        nmeas = flex.double(hklin_df["nmeas"])  # ASK double x int
        m = miller.array(
            miller_set=miller.set(cs, indices),
            data=nmeas,
        )
    else:
        m = None
        # TO DO warning
    m = m.resolution_filter(d_max=d_max, d_min=d_min)
    return m


def calc_stats_merged(m_all_i, m_all_nmeas, d_max=0, d_min=0, n_bins=10):
    stats = {"overall": {}, "binned": {}}
    res_low, res_high = m_all_i.d_max_min()
    n_unique = m_all_i.size()
    completeness = m_all_i.as_non_anomalous_set().completeness()
    m_all_i = m_all_i.map_to_asu()
    m_all_i.setup_binner(n_bins=n_bins)
    m_all_i = m_all_i.sort("packed_indices")
    m_all_i.setup_binner(n_bins=n_bins) # need to redo the binner as is lost above

    # overall values
    i_mean = m_all_i.mean()
    i_sig = m_all_i.i_over_sig_i()
    # n_obs = int(m_all_nmeas.sum())
    n_obs = sum(m_all_nmeas.data().iround())
    multiplicity = m_all_nmeas.mean()
    stats["overall"]["d_max"] = round(res_low, 3)
    stats["overall"]["d_min"] = round(res_high, 3)
    stats["overall"]["n_unique"] = n_unique
    stats["overall"]["n_obs"] = n_obs
    stats["overall"]["completeness"] = round(completeness * 100, 2)
    stats["overall"]["multiplicity"] = round(multiplicity, 2)
    stats["overall"]["I"] = round(i_mean, 2)
    stats["overall"]["IsigI"] = round(i_sig, 2)
    print(f"#observed: {n_obs}")
    print(f"#unique: {m_all_i.size()}")
    print(f"completeness = {completeness * 100:.2f} %")
    print(f"multiplicity = {multiplicity:.2f}")
    print(f"<I> = {i_mean:.1f}")
    print(f"<I/sigma(I)> = {i_sig:.1f}")

    # binned values
    stats["binned"]["d_max"] = []
    stats["binned"]["d_min"] = []
    stats["binned"]["n_obs"] = []
    stats["binned"]["n_unique"] = []
    stats["binned"]["completeness"] = []
    stats["binned"]["multiplicity"] = []
    stats["binned"]["I"] = []
    stats["binned"]["IsigI"] = []
    # print("Binned values:")
    # print(f"d_max   d_min  #unique %complet. multip. <I> <I/sigma(I)>")
    for i_bin in m_all_i.binner().range_used():
        sel = m_all_i.binner().selection(i_bin)
        m_all_i_sel = m_all_i.select(sel)
        m_all_nmeas_sel = m_all_nmeas.select(sel)
        res_low, res_high = m_all_i_sel.d_max_min()
        n_unique = m_all_i_sel.size()
        completeness = m_all_i_sel.as_non_anomalous_set().completeness(d_max=res_low)
        # n_obs = int(m_all_nmeas_sel.sum())
        n_obs = sum(m_all_nmeas_sel.data().iround())
        n_ref_nmeas = int(m_all_nmeas_sel.size())
        multiplicity = m_all_nmeas_sel.mean()
        i_mean =  m_all_i_sel.mean()
        i_sig =  m_all_i_sel.i_over_sig_i()
        stats["binned"]["d_max"].append(round(res_low, 3))
        stats["binned"]["d_min"].append(round(res_high, 3))
        stats["binned"]["n_obs"].append(n_obs)
        stats["binned"]["n_unique"].append(n_unique)
        stats["binned"]["completeness"].append(round(completeness * 100, 2))
        stats["binned"]["multiplicity"].append(round(multiplicity, 2))
        stats["binned"]["I"].append(round(i_mean, 2))
        stats["binned"]["IsigI"].append(round(i_sig, 2))
        # print(f"{res_low:.3f}  {res_high:.3f}  {n_unique}   {completeness * 100:.1f}     {multiplicity:.1f}  {i_mean:.2f}  {i_sig:.2f}")
        if n_unique != n_ref_nmeas:
            sys.stderr.write(
                f"WARNING: Numbers of reflections in bins for the intensity Miller"
                f"array and multiplicity Miller array do not equal.\n"
                f"{n_ref} is not {n_ref_nmeas} in bin {round(res_low, 3)}-{round(res_high, 3)}\n")
    return stats


def calc_stats_compare(m_all, m1, m2, d_max, d_min, n_bins):
    stats = {"overall": {}, "binned": {}}
    m1 = m1.resolution_filter(d_max=d_max, d_min=d_min)
    m1 = m1.map_to_asu()
    m1 = m1.sort("packed_indices")

    m2 = m2.resolution_filter(d_max=d_max, d_min=d_min)
    m2 = m2.map_to_asu()
    m2 = m2.sort("packed_indices")

    m1, m2 = m1.common_sets(m2, assert_no_singles=False)

    # overall values   
    cc = flex.linear_correlation(m1.data(), m2.data()).coefficient()
    CCstar = calc_CCstar(cc)
    rsplit = calc_rsplit(m1.data(), m2.data())
    stats["overall"]["cc"] = round(cc, 3)
    stats["overall"]["CCstar"] = round(CCstar, 3)
    stats["overall"]["rsplit"] = round(rsplit, 3)
    # print("Overall values:")
    print(f"CC1/2 = {cc:.3f}\nCC* = {CCstar:.3f}\nRsplit = {rsplit:.3f}")

    # binned values
    m1.use_binning(m_all.binner())
    m2.use_binning(m_all.binner())
    stats["binned"]["cc"] = []
    stats["binned"]["CCstar"] = []
    stats["binned"]["rsplit"] = []
    # print("Binned values:")
    # print(f"d_max   d_min  CC1/2  CC*   Rsplit")
    for i_bin in m1.binner().range_used():
        sel = m1.binner().selection(i_bin)
        m1_sel = m1.select(sel)
        m2_sel = m2.select(sel)
        res_low, res_high = m1_sel.d_max_min()
        cc = flex.linear_correlation(m1_sel.data(), m2_sel.data()).coefficient()
        CCstar = calc_CCstar(cc)
        rsplit = calc_rsplit(m1_sel.data(), m2_sel.data())
        stats["binned"]["cc"].append(round(cc, 3))
        stats["binned"]["CCstar"].append(round(CCstar, 3))
        stats["binned"]["rsplit"].append(round(rsplit, 3))
        # print(f"{res_low:.3f}  {res_high:.3f}  {cc:.3f} {CCstar:.3f} {rsplit:.3f}")
    return stats


def calc_cc_rsplit(half_dataset, spacegroup, cell, d_max=0, d_min=0, n_bins=10):
    """Code from James"""
    from cctbx import miller, crystal, uctbx, sgtbx, xray
    from cctbx.array_family import flex
    import pandas as pd
    import math
    
    # read the data from the text file to get lists of the values
    # this code assumes the miller indices are the same and in the same order between files
    # TO DO: I/sigma?
    # expected format of a fixed-width .hkl file:
    # it contains just data, no header or footer
    #    h    k    l          I    phase   sigma(I)   nmeas
    half1 = pd.read_csv(
        half_dataset[0], header=None, index_col=False, sep=r'\s+',
        names=("h", "k", "l", "I", "phase", "sigma(I)", "nmeas"))
    half2 = pd.read_csv(
        half_dataset[1], header=None, index_col=False, sep=r'\s+',
        names=("h", "k", "l", "I", "phase", "sigma(I)", "nmeas"))
    h1 = flex.int(half1["h"])
    k1 = flex.int(half1["k"])
    l1 = flex.int(half1["l"])
    I1 = flex.double(half1["I"])
    sig1 = flex.double(half1["sigma(I)"])
    nmeas1 = flex.double(half1["nmeas"])
    h2 = flex.int(half2["h"])
    k2 = flex.int(half2["k"])
    l2 = flex.int(half2["l"])
    I2 = flex.double(half2["I"])
    sig2 = flex.double(half2["sigma(I)"])

    cs = crystal.symmetry(
        unit_cell=uctbx.unit_cell(cell),
        space_group=sgtbx.space_group_info(spacegroup).group()
    )

    indices1 = flex.miller_index(h1,k1,l1)
    m1 = miller.array(
        miller_set=miller.set(cs, indices1),
        data=I1,
        sigmas=sig1,
    )
    m1_nmeas = miller.array(
        miller_set=miller.set(cs, indices1),
        data=nmeas1,
    )
    m1 = m1.resolution_filter(d_max=d_max, d_min=d_min)
    m1.set_observation_type_xray_intensity()
    m1 = m1.map_to_asu()
    m1 = m1.sort("packed_indices")
    m1_nmeas = m1_nmeas.resolution_filter(d_max=d_max, d_min=d_min)
    #m1_nmeas = m1_nmeas.map_to_asu()
    #m1_nmeas = m1_nmeas.sort("packed_indices")
    
    indices2 = flex.miller_index(h2,k2,l2)
    m2 = miller.array(
        miller_set=miller.set(cs, indices2),
        data=I2,
        sigmas=sig2,
    )
    m2 = m2.resolution_filter(d_max=d_max, d_min=d_min)
    m1.set_observation_type_xray_intensity()
    m2 = m2.map_to_asu()
    m2 = m2.sort("packed_indices")
    m1, m2 = m1.common_sets(m2, assert_no_singles=False)
    #m1_nmeas, trash = m1_nmeas.common_sets(m1, assert_no_singles=False)
    print(f"{m1.size()} {m1_nmeas.size()}")

    # overall values   
    cc = flex.linear_correlation(m1.data(), m2.data()).coefficient()
    rsplit = calc_rsplit(m1.data(), m2.data())
    i_sig = m1.i_over_sig_i()
    multiplicity = m1_nmeas.mean()
    print("\nOverall values:")
    print(f"CC1/2 = {cc:.3f} \nRsplit = {rsplit:.3f} \nI/sigma(I) = {i_sig:.3f}")
    print(f"Multiplicity = {multiplicity:.1f}")
    print(f"{m1.size()} {m1_nmeas.size()}")
    print("")

    # binned values
    m1.setup_binner(n_bins=n_bins)
    m2.use_binning(m1.binner())
    print("\nBinned values:")
    print(f"d_max   d_min   #refl CC1/2  Rsplit") #"  I/sigma(I)")
    for i_bin in m1.binner().range_used():
        sel = m1.binner().selection(i_bin)
        m1_sel = m1.select(sel)
        #m1_sel.show_array()
        m2_sel = m2.select(sel)
        res_low, res_high = m1_sel.d_max_min()
        n_ref = m1_sel.size()
        cc = flex.linear_correlation(m1_sel.data(), m2_sel.data()).coefficient()
        rsplit = calc_rsplit(m1_sel.data(), m2_sel.data())
        i_sig = m1_sel.i_over_sig_i()
        print(f"{res_low:.3f}  {res_high:.3f}  {n_ref}  {cc:.3f}  {rsplit:.3f}")#  {i_sig:.3f}")


def stats_to_xml(stats):  #, xmlout="program.xml"):
    lines = []
    lines.append("<import_serial>")
    for key1, key2 in stats.items():
        lines.append(f"\t<{key1}>")  # overall or binned
        if key1 == "overall":
            for key_2, value in key2.items():
                lines.append(f"\t\t<{key_2}>{value}</{key_2}>")
        elif key1 == "binned":
            for i in range(len(list(key2.values())[0])):  # for individual bins
                lines.append(f"\t\t<bin>")
                lines.append(f"\t\t\t<n_bin>{i + 1}</n_bin>")
                for key_2, key_3 in key2.items():
                    lines.append(f"\t\t\t<{key_2}>{key_3[i]}</{key_2}>")  # statistic
                    if key_2 == "d_min":
                        over_d_min_sq = 1 / (key_3[i] * key_3[i])
                        over_d_min_sq = round(over_d_min_sq, 4)
                        lines.append(f"\t\t\t<one_over_d_min_sq>{over_d_min_sq}</one_over_d_min_sq>")
                lines.append(f"\t\t</bin>")
            #for key_2, key_3 in key2.items():
            #    lines.append(f"\t\t<{key_2}>")  # statistic
            #    for i, value in enumerate(key_3):  # key_3 is list
            #        lines.append(f"\t\t\t<bin" + str(i+1) + f">{value}</bin" + str(i+1) + ">")
            #    lines.append(f"\t\t</{key_2}>")
        lines.append(f"\t</{key1}>")
    lines.append("</import_serial>")
    stats_xml = '\n'.join(lines)
    return stats_xml


def stats_binned_print(stats_binned):  # , half_dataset_available=False):
    # stats_print = "%9s%9s%11s%8s%9s%9s%9s%9s"
    header = ["d_max", "d_min", "#obs", "#uniq", "mult.", "%comp", "<I>", "<I/sI>"]
    stats_print_format_header = '{:>8}{:>8}{:>9}{:>8}{:>8}{:>8}{:>9}{:>8}'
    stats_print_format_values =  '{:>8.2f}{:>8.2f}{:>9d}{:>8d}{:>8.2f}{:>8.2f}{:>9.1f}{:>8.1f}'
    # if half_dataset_available:
    if len(stats_binned) == 11:
        # stats_print += "%9s%9s%9s"
        header += ["cc1/2", "cc*", "r_split"]
        stats_print_format_header += '{:>8}{:>8}{:>8}'
        stats_print_format_values += '{:>8.3f}{:>8.3f}{:>8.3f}'
    # print(stats_print % tuple(stats_print_header))
    print(stats_print_format_header.format(*header))
    for i in range(len(stats_binned["d_max"])):
        values = []
        values.append(stats_binned["d_max"][i])
        values.append(stats_binned["d_min"][i])
        values.append(stats_binned["n_obs"][i])
        values.append(stats_binned["n_unique"][i])
        values.append(stats_binned["multiplicity"][i])
        values.append(stats_binned["completeness"][i])
        values.append(stats_binned["I"][i])
        values.append(stats_binned["IsigI"][i])
        # if half_dataset_available:
        if len(stats_binned) == 11:
            values.append(stats_binned["cc"][i])
            values.append(stats_binned["CCstar"][i])
            values.append(stats_binned["rsplit"][i])
        # print(stats_print % tuple(values))
        print(stats_print_format_values.format(*values))
    return


class MyArgumentParser(argparse.ArgumentParser):
    """ Helper class for `argparse`

    It adds an attribute `add_argument_with_check` that check if file
    given in argument exist but not open it.
    Inspired by `https://codereview.stackexchange.com/questions/28608/
    checking-if-cli-arguments-are-valid-files-directories-in-python`
    """
    def __is_valid_file(self, parser, arg):
        """Checks if file
        given in argument `arg` exists but does not open it.
        If not, abort.

        Args:
            self
            parser: parser of `argparse`
            arg (str): argument of `argparse`

        Returns:
            str: Name of the checked file
        """
        if not os.path.isfile(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            # File exists so return the filename
            return arg

    def add_argument_with_check(self, *args, **kwargs):
        """New attribute for `argparse` that checks if file
        given in argument exist but does not open it.
        """
        # Look for your FILE settings
        # type = lambda x: self.__is_valid_file(self, x) # PEP8 E731
        def type(x):
            return self.__is_valid_file(self, x)
        kwargs['type'] = type
        self.add_argument(*args, **kwargs)

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        print('Run `ccp4-python -m import_serial -h` to show the help message.')
        sys.exit(2)


def get_cell_cellfile(cellfile):
    with open(cellfile, 'r') as f:
        lines = f.readlines()
    cell = [None, None, None, None, None, None]
    cell_string = None
    for line in lines:
        try:
            if line.split()[0] == "al" and line.split()[-1] == "deg":
                cell[3] = float(line.split()[-2])
            elif line.split()[0] == "be" and line.split()[-1] == "deg":
                cell[4] = float(line.split()[-2])
            elif line.split()[0] == "ga" and line.split()[-1] == "deg":
                cell[5] = float(line.split()[-2])
            elif line.split()[0] == "a" and line.split()[-1] == "A":
                cell[0] = float(line.split()[-2])
            elif line.split()[0] == "b" and line.split()[-1] == "A":
                cell[1] = float(line.split()[-2])
            elif line.split()[0] == "c" and line.split()[-1] == "A" and not "centering" in line:
                cell[2] = float(line.split()[-2])
        except:
            continue
    print("")
    if (cell[0] and cell[1] and cell[2] and cell[3] and cell[4] and cell[5]):
        cell_string = " ".join(map(str, cell))
        print(f"Unit cell parameters found in file {cellfile}:")
        print(cell_string)
    else:
        sys.stderr.write(
            f"WARNING: Unit cell parameters could not be parsed from "
            f"the file {cellfile}.\n"
            f"Attempt to find the unit cell parameters found "
            f"in this file: " + " ".join(map(str, cell)) + "\n")
        cell = None
    return cell, cell_string


def get_cell_streamfile(streamfile):
    cell = [None, None, None, None, None, None]
    cell_string = None
    if os.path.isfile(streamfile + "_tmp"):
        os.remove(streamfile + "_tmp")
    with open(streamfile, "r") as file1:
        with open(streamfile + "_tmp", "a+") as file2:
            for line in file1:
                if "Cell parameters " in line and len(line.split()) == 10:
                    file2.write(line)
    if os.path.isfile(streamfile + "_tmp"):
        cell_df = pd.read_csv(
            streamfile + "_tmp", header=None, index_col=False, sep=r'\s+',
            names=("none1", "none2", "a", "b", "c", "none3", "alpha", "beta", "gamma", "none4"))
        cell_df = cell_df.drop(columns=["none1", "none2", "none3", "none4"])
        cell_df = cell_df.astype(float)
        cell[0] = round(cell_df.mean()["a"] * 10, 2)
        cell[1] = round(cell_df.mean()["b"] * 10, 2)
        cell[2] = round(cell_df.mean()["c"] * 10, 2)
        cell[3] = round(cell_df.mean()["alpha"], 2)
        cell[4] = round(cell_df.mean()["beta"], 2)
        cell[5] = round(cell_df.mean()["gamma"], 2)
        cell_string = " ".join(map(str, cell))
        print("")
        print(f"Unit cell parameters fit using file {streamfile}:")
        print(cell_string)
        os.remove(streamfile + "_tmp")
    else:
        sys.stderr.write(
            f"WARNING: Unit cell parameters could not be fitted from "
            f"the file {streamfile}.\n")
        cell = None
    return cell, cell_string


def get_wavelength_streamfile(streamfile):
    energy_eV = None
    wavelength = None
    if os.path.isfile(streamfile + "_tmp"):
        os.remove(streamfile + "_tmp")
    with open(streamfile, "r") as file1:
        with open(streamfile + "_tmp", "a+") as file2:
            for line in file1:
                if "photon_energy_eV" in line and len(line.split()) == 3:
                    file2.write(line)
    if os.path.isfile(streamfile + "_tmp"):
        energy_eV_df = pd.read_csv(
            streamfile + "_tmp", header=None, index_col=False, sep=r'\s+',
            names=("none1", "none2", "energy_eV"))
        energy_eV_df = energy_eV_df.drop(columns=["none1", "none2"])
        energy_eV_df = energy_eV_df.astype(float)
        energy_eV = energy_eV_df.median()["energy_eV"]
        wavelength = 12398.425 / energy_eV
        wavelength = round(wavelength, 5)
        print("")
        print(f"Wavelength median using file {streamfile}:")
        print(str(wavelength))
        os.remove(streamfile + "_tmp")
    else:
        sys.stderr.write(
            f"WARNING: Wavelength could not be fitted from "
            f"the file {streamfile}.\n")
        wavelength = 0
    return wavelength


def get_wavelength_reference(ref):
    if reflection_file_reader.any_reflection_file(ref).file_type() == 'ccp4_mtz':
        try:
            mtz_object = mtz.object(file_name=ref)
            crystal = mtz.crystal(mtz_object=mtz_object, i_crystal=1)
            dataset = mtz.dataset(mtz_crystal=crystal, i_dataset=0)
            wavelength = round(float(dataset.wavelength()), 5)
            print("")
            print(f"Wavelength found in {ref}:")
            print(str(wavelength))
        except:
            sys.stderr.write(
                f"WARNING: Wavelength could not found in "
                f"the file {ref}.\n")
            wavelength = 0
    else:
        wavelength = 0
    return wavelength


def get_cs_reference(reference):
    from iotbx import file_reader
    file = file_reader.any_file(reference)
    try:
        cs = file.crystal_symmetry()
        spacegroup = file.crystal_symmetry().space_group().info()
        cell = file.crystal_symmetry().unit_cell().parameters()
        cell_string = " ".join(map(str, cell))
        print("")
        print(f"Symmetry from the reference file {reference}:")
        print(str(cs))
    except NotImplementedError:
        sys.stderr.write(
            f"WARNING: Symmetry could not be found in the provided "
            f"reference file {reference}.\n")
        return None, None, None
    return cs, spacegroup, cell_string


def run():
    from . import __version__
    # if not which("f2mtz"):
    #     sys.stderr.write(f"ERROR: Program f2mtz from CCP4 is not available.\n"
    #                      "Did you source the paths to CCP4 executables?"
    #                      "Aborting.\n")
    #     sys.exit(1)
    parser = MyArgumentParser(
        description="Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4"
    )
    parser.add_argument_with_check(
        "--hklin", "--HKLIN",
        help="Specify merged mtz file from xia2.ssx or merged hkl file from CrystFEL",
        type=str,
        required=True
    )
    parser.add_argument(
        "--spacegroup",
        type=str,
        help="Space group",
    )
    parser.add_argument(
        "--cell",
        type=float,
        nargs=6,
        help="Unit cell parameters divided by spaces, e.g. 60 50 40 90 90 90",
        metavar=("a", "b", "c", "alpha", "beta", "gamma"),
    )
    parser.add_argument_with_check(
        "--cellfile",
        help="Cell file from CrystFEL",
        type=str,
    )
    parser.add_argument_with_check(
        "--streamfile",
        help="Stream file from CrystFEL",
        type=str,
    )
    parser.add_argument_with_check(
        "--reference", "--ref", "--pdb", "--cif", "--mmcif",
        metavar="REFERENCE",
        help="Reference file (PDB, mmCIF or MTZ) to provide spacegroup and unit cell",
        type=str,
        dest="ref",
    )
    parser.add_argument_with_check(
        "--half-dataset",
        metavar=("HKL1", "HKL2"),
        help="CrystFEL only: two half-data-set merge files (usually .hkl1 and .hkl2)",
        type=str,
        nargs=2,
    )
    parser.add_argument(
        "--project",
        type=str,
        help="Project name",
    )
    parser.add_argument(
        "--crystal",
        type=str,
        help="Crystal name",
        dest='cryst',
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Dataset name",
    )
    parser.add_argument(
        "--dmax", "--lowres",
        type=float,
        help="Low-resolution cutoff",
        default=0,
        dest='d_max',
    )
    parser.add_argument(
        "--dmin", "--highres",
        type=float,
        help="High-resolution cutoff",
        default=0,
        dest='d_min',
    )
    parser.add_argument(
        "--nbins", "--nshells",
        type=int,
        help="Number of resolution bins",
        default=10,
        dest='n_bins',
    )
    parser.add_argument(
        "--wavelength", "-w",
        type=float,
        help="Wavelength (only for data from CrystFEL)",
        default=0,
    )
    # TO DO: check whether the files exist etc.
    args = parser.parse_args()

    print("")
    print("Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4")
    print("")
    print("Command line arguments:")
    print(" ".join(sys.argv[1:]))
    print("")
    print("Input parameters:")
    for arg in vars(args):
        if getattr(args, arg):
            print('  {} {}'.format(arg, getattr(args, arg) or ''))

    hklin = args.hklin
    if reflection_file_reader.any_reflection_file(hklin).file_type() == 'ccp4_mtz':
        hklin_format = "dials"
    elif reflection_file_reader.any_reflection_file(hklin).file_type() == None:
        hklin_format = "crystfel"
        spacegroup = None
        cell = None

    if not args.project:
        project = "project"
    else:
        project = args.project
    if not args.cryst:
        cryst = "crystal"
    else:
        dataset = args.dataset
    if not args.dataset:
        dataset = "dataset"
    else:
        dataset = args.dataset
    prefix = f"{project}_{dataset}"
    prefix = "".join(i for i in prefix if i not in r"\/:*?<>|")
    hklout = f"{prefix}.mtz"
    logout = f"{prefix}.html"
    jsonout = f"{prefix}.json"
    xmlout = "program.xml"
    # xmlout = f"{prefix}.xml"
    d_max = args.d_max
    d_min = args.d_min
    n_bins = args.n_bins
    wavelength = args.wavelength

    # process symmetry: space group and unit cell parameters
    cs = None
    if args.spacegroup:
        if args.cell:
            cell = args.cell
        elif args.cellfile:
            cell, cell_string = get_cell_cellfile(args.cellfile)
        elif args.streamfile:
            cell, cell_string = get_cell_streamfile(args.streamfile)
        if args.cell or args.cellfile or args.streamfile:  # everything except reference file
            spacegroup = args.spacegroup
            cs = crystal.symmetry(
                unit_cell=uctbx.unit_cell(cell),
                space_group=sgtbx.space_group_info(spacegroup).group())
        elif args.ref:
            cs, spacegroup, cell_string = get_cs_reference(args.ref)
    elif args.ref:
        cs, spacegroup, cell_string = get_cs_reference(args.ref)
    # crystfel: check whether we know spacegroup and cell
    if hklin_format == "crystfel" and not cs:
        # raise error and abort
        if not args.spacegroup and (args.cell or args.cellfile):
            # error missing spacegroup
            sys.stderr.write(
                "ERROR: Space group was not specified but is required for CrystFEL.\n"
                "Specify space group explicitly (option --spacegroup) "
                "or provide reference PDB, mmCIF or MTZ file (option --reference).\n")
        elif (not args.cell or not args.cellfile) and args.spacegroup:
            # error missing cell
            sys.stderr.write(
                "ERROR: Unit cell parameters were not specified but are required for CrystFEL.\n"
                "Specify unit cell parameters explicitly (options  --cell or --cellfile) "
                "or provide reference PDB, mmCIF or MTZ file (option --reference).\n")
        else:  # if not args.spacegroup and not args.cell and not args.cellfile and not args.ref:
            # error missing everything
            sys.stderr.write(
                "ERROR: Unit cell parameters and spacegroup were not specified but are required for CrystFEL.\n"
                "Specify unit cell parameters (options  --cell or --cellfile) "
                "and space group (option --spacegroup) "
                "or provide reference PDB, mmCIF or MTZ file (option --reference).\n")
        sys.stderr.write("Aborting.\n")
        sys.exit(1)
    if hklin_format == "crystfel" and args.streamfile and not wavelength:
        wavelength = get_wavelength_streamfile(args.streamfile)
    elif hklin_format == "crystfel" and args.ref and not wavelength:
        wavelength = get_wavelength_reference(args.ref)

    print("")
    print("")
    print("DATA STATISTICS:")
    print("================")
    print("")
    try:
        m1 = None
        m2 = None
        # load data to Miller arrays
        if hklin_format == "dials":
            # miller_arrays = reflection_file_reader.any_reflection_file(file_name=hklin).as_miller_arrays()
            miller_arrays = mtz.object(hklin).as_miller_arrays()
            m_all_i = None
            m1_all_nmeas = None
            m2_all_nmeas = None
            for i, column in enumerate(miller_arrays):
                # print(str(column.info().labels))
                # if column.is_xray_intensity_array() and str(column.info()).split(".mtz:")[1].split(",")[0] == "IMEAN":
                # elif column.is_xray_intensity_array() and str(column.info()).split(".mtz:")[1].split(",")[0] == "IHALF1":
                # elif column.is_xray_intensity_array() and str(column.info()).split(".mtz:")[1].split(",")[0] == "IHALF2":
                if column.is_xray_intensity_array() and column.info().labels == ["IMEAN", "SIGIMEAN"]:
                    m_all_i = column
                elif column.is_xray_intensity_array() and column.info().labels == ["IHALF1", "SIGIHALF1"]:
                    m1 = column
                elif column.is_xray_intensity_array() and column.info().labels == ["IHALF2", "SIGIHALF2"]:
                    m2 = column
                elif column.info().labels == ['N']:
                    m_all_nmeas = column.as_double()
        elif hklin_format == "crystfel":
            if args.half_dataset:
                half_dataset = args.half_dataset
            elif os.path.isfile(hklin) and os.path.isfile(hklin + "1") and os.path.isfile(hklin + "2"):
                half_dataset = (hklin + "1", hklin + "2")
                print(
                    f"Half-dataset files were found automatically and will be used "
                    f"for calculation of statistics: {half_dataset[0]} {half_dataset[1]}")
            else:
                half_dataset = None
            m_all_i = get_miller_array_crystfel(hklin, cs, "I", d_max=d_max, d_min=d_min)
            m_all_nmeas = get_miller_array_crystfel(hklin, cs, "nmeas", d_max=d_max, d_min=d_min)
            if half_dataset:
                m1 = get_miller_array_crystfel(half_dataset[0], cs, "I", d_max=d_max, d_min=d_min)
                m2 = get_miller_array_crystfel(half_dataset[1], cs, "I", d_max=d_max, d_min=d_min)

        # set d_min, d_max and binning to miller arrays
        m_all_i = m_all_i.resolution_filter(d_max=d_max, d_min=d_min)
        m_all_nmeas = m_all_nmeas.resolution_filter(d_max=d_max, d_min=d_min)
        m_all_i.setup_binner(n_bins=n_bins)
        m_all_nmeas.use_binning(m_all_i.binner())

        # calculate and print statistics
        print("Overall values:\n")
        stats_merged = calc_stats_merged(m_all_i, m_all_nmeas, d_max, d_min, n_bins)
        if m1 and m2:
            # calculate statistics CC1/2, CC* and Rsplit
            stats_compare = calc_stats_compare(m_all_i, m1, m2, d_max, d_min, n_bins)
            stats_overall = {**stats_merged["overall"], **stats_compare["overall"]}
            stats_binned = {**stats_merged["binned"], **stats_compare["binned"]}
        else:
            stats_overall = {**stats_merged["overall"]}
            stats_binned = {**stats_merged["binned"]}
        print("\nBinned values:\n")
        stats_binned_print(stats_binned)
        
        # save statistics to files
        stats = {"overall": stats_overall, "binned": stats_binned}
        stats_json = json.dumps(stats, indent=4)
        stats_xml = stats_to_xml(stats)  #, xmlout)
        with open(jsonout, "w") as f:
            f.write(stats_json)
        with open(xmlout, "w") as f:
            f.write(stats_xml)
    except RuntimeError:
        traceback.print_exc()
        sys.stderr.write("WARNING: Statistics could not be calculated.\n")
    # TO DO check if output files with same name exists already

    # compare_hkl $inp1 $inp2 -y $symm -p $pdb --fom=$mode --highres=$highres --nshells=20 --shell-file="stat/${basename}-$mode".dat 2>>stat/${basename}.log
    # check_hkl -p $pdb --nshells=20 --highres=$highres -y $pg  --shell-file="stat/${basename}-shells".dat $inp 2>>stat/${basename}.log

    #hkltmp = hkl_strip(args.hklin)
    #command = ["f2mtz", "HKLIN", hkltmp, "HKLOUT", hklout]
    #inp = f"""TITLE Reflections from CrystFEL
    #NAME PROJECT {project} CRYSTAL {cryst} DATASET {dataset}
    #CELL {cell_str}
    #SYMM {spacegroup}
    #LABOUT H K L IMEAN SIGIMEAN
    #CTYPE  H H H J     Q
    #FORMAT '(3(F4.0,1X),F10.2,10X,F10.2)'
    #END
    #"""
    #with open(logout, "w") as logfile:
    #    p = subprocess.Popen(
    #        command, stdin=subprocess.PIPE,
    #        stdout=logfile, stderr=logfile,
    #        encoding="utf-8")  # shell=settings["sh"])
    #    output, err = p.communicate(inp)
    #    # rc = p.returncode
    # remove hkltmp ?
    # print(f"MTZ file created: {hklout}")
    if hklin_format == "crystfel":
        mtz_dataset = m_all_i.as_mtz_dataset(column_root_label="IMEAN", wavelength=wavelength)
        mtz_dataset.add_miller_array(m_all_nmeas, column_root_label="NMEAS")
        # mtz_dataset.add_miller_array(r_free_flags, column_root_label="FreeR_flag")
        mtz_dataset.mtz_object().write(file_name=hklout)
        print(f"\nMTZ file created: {hklout}")
        # Clean up a bit
        files_to_remove = [hklin + "_tmp"]
        if half_dataset:
            files_to_remove += [f'{half_dataset[0]}_tmp']
            files_to_remove += [f'{half_dataset[1]}_tmp']
        for f in files_to_remove:
            if os.path.isfile(f):
                os.remove(f)
    elif hklin_format == "dials":
        import shutil
        shutil.copy2(hklin, hklout)
    return
