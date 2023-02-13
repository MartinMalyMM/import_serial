# author: Martin Maly
# email: martin.maly.mm@email.cz
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
from cctbx import miller, crystal, uctbx, sgtbx, xray
from cctbx.array_family import flex
# TO DO: 3 input files
# TO DO: unit cell parameters from stream file
# TO DO: check output MTZ
# TO DO: wavelength


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


def get_miller_array(hklin, cs, values="I", d_max=0, d_min=0):
    assert values == "I" or values == "nmeas"
    # read data from the text file
    # expected format of a fixed-width .hkl file:
    #    h    k    l          I    phase   sigma(I)   nmeas
    hklin_strip = hkl_strip(hklin)
    hklin_df = pd.read_csv(
        hklin_strip, header=None, index_col=False, sep='\s+',
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
    m_all_i = m_all_i.resolution_filter(d_max=d_max, d_min=d_min)
    res_low, res_high = m_all_i.d_max_min()
    n_ref = m_all_i.size()
    completeness = m_all_i.as_non_anomalous_set().completeness()
    m_all_i = m_all_i.map_to_asu()
    m_all_i = m_all_i.sort("packed_indices")
    m_all_nmeas = m_all_nmeas.resolution_filter(d_max=d_max, d_min=d_min)

    # overall values
    i_mean = m_all_i.mean()
    i_sig = m_all_i.i_over_sig_i()
    multiplicity = m_all_nmeas.mean()
    stats["overall"]["d_max"] = round(res_low, 3)
    stats["overall"]["d_min"] = round(res_high, 3)
    stats["overall"]["n_unique"] = n_ref
    stats["overall"]["completeness"] = round(completeness * 100, 1)
    stats["overall"]["multiplicity"] = round(multiplicity, 1)
    stats["overall"]["I"] = round(i_mean, 2)
    stats["overall"]["IsigI"] = round(i_sig, 2)
    print("Overall values:")
    print(f"#unique: {m_all_i.size()}")
    print(f"completeness = {completeness * 100:.1f} %")
    print(f"multiplicity = {multiplicity:.1f}")
    print(f"<I> = {i_mean:.2f}")
    print(f"<I/sigma(I)> = {i_sig:.2f}")
    print("")

    # binned values
    m_all_i.setup_binner(n_bins=n_bins)  # TO DO when/where to set up binner?
    m_all_nmeas.use_binning(m_all_i.binner())
    stats["binned"]["d_max"] = []
    stats["binned"]["d_min"] = []
    stats["binned"]["n_unique"] = []
    stats["binned"]["completeness"] = []
    stats["binned"]["multiplicity"] = []
    stats["binned"]["I"] = []
    stats["binned"]["IsigI"] = []
    print("Binned values:")
    print(f"d_max   d_min  #unique %complet. multip. <I> <I/sigma(I)>")
    for i_bin in m_all_i.binner().range_used():
        sel = m_all_i.binner().selection(i_bin)
        m_all_i_sel = m_all_i.select(sel)
        m_all_nmeas_sel = m_all_nmeas.select(sel)
        res_low, res_high = m_all_i_sel.d_max_min()
        n_ref = m_all_i_sel.size()
        completeness = m_all_i_sel.as_non_anomalous_set().completeness(d_max=res_low)
        n_ref_nmeas = m_all_nmeas_sel.size()
        multiplicity = m_all_nmeas_sel.mean()
        i_mean =  m_all_i_sel.mean()
        i_sig =  m_all_i_sel.i_over_sig_i()
        stats["binned"]["d_max"].append(round(res_low, 3))
        stats["binned"]["d_min"].append(round(res_high, 3))
        stats["binned"]["n_unique"].append(n_ref)
        stats["binned"]["completeness"].append(round(completeness * 100, 1))
        stats["binned"]["multiplicity"].append(round(multiplicity, 1))
        stats["binned"]["I"].append(round(i_mean, 2))
        stats["binned"]["IsigI"].append(round(i_sig, 2))
        print(f"{res_low:.3f}  {res_high:.3f}  {n_ref}   {completeness * 100:.1f}     {multiplicity:.1f}  {i_mean:.2f}  {i_sig:.2f}")
    if n_ref != n_ref_nmeas:
        sys.stderr.write(
            "WARNING: Numbers of reflections in bins for the intensity Miller"
            "array and multiplicity Miller array do not equal.\n"
            "{n_ref} is not {n_ref_nmeas}\n")
    print("")
    print("")
    return stats


def calc_stats_compare(m1, m2, d_max, d_min, n_bins):
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
    print("Overall values:")
    print(f"CC1/2 = {cc:.3f} \nCC* = {CCstar:.3f}\nRsplit = {rsplit:.3f}")
    print("")

    # binned values
    m1.setup_binner(n_bins=n_bins)
    m2.use_binning(m1.binner())
    stats["binned"]["cc"] = []
    stats["binned"]["CCstar"] = []
    stats["binned"]["rsplit"] = []
    print("Binned values:")
    print(f"d_max   d_min  CC1/2  CC*   Rsplit")
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
        print(f"{res_low:.3f}  {res_high:.3f}  {cc:.3f} {CCstar:.3f} {rsplit:.3f}")
    print("")
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
        half_dataset[0], header=None, index_col=False, sep='\s+',
        names=("h", "k", "l", "I", "phase", "sigma(I)", "nmeas"))
    half2 = pd.read_csv(
        half_dataset[1], header=None, index_col=False, sep='\s+',
        names=("h", "k", "l", "I", "phase", "sigma(I)", "nmeas"))
    #print(str(half1["nmeas"][:20])) ERASE
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
    #multiplicity = m1_nmeas.data().mean() NO
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


def run():
    from . import __version__
    if not which("f2mtz"):
        sys.stderr.write(f"ERROR: Program f2mtz from CCP4 is not available.\n"
                         "Did you source the paths to CCP4 executables?"
                         "Aborting.\n")
        sys.exit(1)
    parser = argparse.ArgumentParser(  # TODO requirement
        description="Covert hkl file(s) from CrystFEL to mtz and calculate statistics"
    )
    parser.add_argument(
        "--hklin", "--HKLIN",
        metavar=("hklin_file"),
        help="Specify one merged hkl file from CrystFEL",
        type=str,
        required=True
    )
    parser.add_argument(
        "--half-dataset",
        metavar=("hkl1", "hkl2"),
        help="Two half-data-set merge files from CrystFEL (usually .hkl1 and .hkl2)",
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
        help="crystal",
        dest='cryst',
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="dataset",
    )
    parser.add_argument(
        "--spacegroup",
        type=str,
        help="Specify space group",
        metavar="spacegroup",
        required=True
    )
    parser.add_argument(
        "--cell",
        type=float,
        nargs=6,
        help="Specify unit cell parameters divided by spaces, e.g. 60 50 40 90 90 90",
        metavar=("cell_a", "cell_b", "cell_c", "cell_alpha", "cell_beta", "cell_gamma"),
        required=True
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
    # TO DO: find cell and symmetry from a reference PDB file
    # TO DO: check whether the files exist etc.
    args = parser.parse_args()

    print("Conversion of CrystFEL HKL file to MTZ and import into CCP4")
    print("")
    print("Command line arguments:")
    print(" ".join(sys.argv[1:]))
    print("")
    print("Input parameters:")
    for arg in vars(args):
        if getattr(args, arg):
            print('  {} {}'.format(arg, getattr(args, arg) or ''))
    # check whether the files exists TO DO
    if args.spacegroup:
        spacegroup = args.spacegroup  # TO DO
    print("")
    if args.cell:
        cell = args.cell  # TO DO

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
    prefix = "".join(i for i in prefix if i not in "\/:*?<>|")
    hklout = f"{prefix}.mtz"
    logout = f"{prefix}.html"
    jsonout = f"{prefix}.json"
    xmlout = "program.xml"
    # xmlout = f"{prefix}.xml"
    d_max = args.d_max
    d_min = args.d_min
    n_bins = args.n_bins
    cell_str = str(args.cell[0]) + " " + str(args.cell[1]) + " " + str(args.cell[2]) + " " + str(args.cell[3]) + " " + str(args.cell[4]) + " " + str(args.cell[5])

    if args.half_dataset:
        half_dataset = args.half_dataset
    elif os.path.isfile(args.hklin) and os.path.isfile(args.hklin + "1") and os.path.isfile(args.hklin + "2"):
        half_dataset = (args.hklin + "1", args.hklin + "2")
        print(
            f"Half-dataset files were found automatically and will be used "
            f"for calculation of statistics: {half_dataset[0]} {half_dataset[1]}")
    # TO DO: check whether the files exists
    # TO DO: parse symmetry from .hkl file if not given
    try:
        print("")
        print("Statistics:")
        print("")
        # load data to Miller array
        cs = crystal.symmetry(
            unit_cell=uctbx.unit_cell(cell),
            space_group=sgtbx.space_group_info(spacegroup).group()
        )
        m_all_i = get_miller_array(args.hklin, cs, "I", d_max=d_max, d_min=d_min)
        m_all_nmeas = get_miller_array(args.hklin, cs, "nmeas", d_max=d_max, d_min=d_min)
        stats_merged = calc_stats_merged(m_all_i, m_all_nmeas, d_max, d_min, n_bins)
        if half_dataset:
            # calculate statistics CC1/2 and Rsplit
            m1 = get_miller_array(half_dataset[0], cs, "I", d_max=d_max, d_min=d_min)
            m2 = get_miller_array(half_dataset[1], cs, "I", d_max=d_max, d_min=d_min)
            stats_compare = calc_stats_compare(m1, m2, d_max, d_min, n_bins)

        stats_overall = {**stats_merged["overall"], **stats_compare["overall"]}
        stats_binned = {**stats_merged["binned"], **stats_compare["binned"]}
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

    mtz_dataset = m_all_i.as_mtz_dataset(column_root_label="IMEAN")
    mtz_dataset.add_miller_array(m_all_nmeas, column_root_label="NMEAS")
    # mtz_dataset.add_miller_array(r_free_flags, column_root_label="FreeR_flag")
    mtz_dataset.mtz_object().write(hklout)
    print(f"MTZ file created: {hklout}")

    # TO DO - analysis as for 'import_merged' - pointless, aimless, phaser_analysis, ctruncate
    return
