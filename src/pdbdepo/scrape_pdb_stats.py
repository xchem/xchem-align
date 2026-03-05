import argparse
import re
from collections import OrderedDict
from pathlib import Path
import gemmi
from gemmi import cif

from pdbdepo import pdb_deposition
from xchemalign import utils

LOG = utils.Logger()


def info(*args, **kwargs):
    if LOG:
        LOG.info(*args, **kwargs)


def warn(*args, **kwargs):
    if LOG:
        LOG.warn(*args, **kwargs)


def error(*args, **kwargs):
    if LOG:
        LOG.error(*args, **kwargs)


refine_ls_shell_regexes = {
    'pdbx_total_number_of_bins_used': r'REMARK\s+\w+\s+TOTAL NUMBER OF BINS USED\s+:\s+([\.\d]+)',
    'd_res_high': r'REMARK\s+\w+\s+BIN RESOLUTION RANGE HIGH\s+:\s+([\.\d]+)',
    'd_res_low': r'REMARK\s+\w+\s+BIN RESOLUTION RANGE LOW\s+:\s+([\.\d]+)',
    'percent_reflns_obs': r'REMARK\s+\w+\s+BIN COMPLETENESS \(WORKING\+TEST\) \(%\)\s+:\s+([\.\d]+)',
    'R_factor_R_work': r'REMARK\s+\w+\s+BIN R VALUE\s+\(WORKING SET\)\s+:\s+([\.\d]+)',
    'R_factor_R_free': r'REMARK\s+\w+\s+BIN FREE R VALUE\s+:\s+([\.\d]+)',
    'number_reflns_R_free': r'REMARK\s+\w+\s+BIN FREE R VALUE SET COUNT\s+:\s+([\.\d]+)',
}

refine_ls_restr_regexes = {
    'r_bond_refined_d': r'REMARK\s+\w+\s+BOND LENGTHS REFINED ATOMS\s+\(A\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_bond_other_d': r'REMARK\s+\w+\s+BOND LENGTHS OTHERS\s+\(A\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_angle_refined_deg': r'REMARK\s+\w+\s+BOND ANGLES REFINED ATOMS\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_angle_other_deg': r'REMARK\s+\w+\s+BOND ANGLES OTHERS\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_dihedral_angle_1_deg': r'REMARK\s+\w+\s+TORSION ANGLES, PERIOD 1\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_dihedral_angle_2_deg': r'REMARK\s+\w+\s+TORSION ANGLES, PERIOD 2\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_dihedral_angle_3_deg': r'REMARK\s+\w+\s+TORSION ANGLES, PERIOD 3\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_dihedral_angle_4_deg': r'REMARK\s+\w+\s+TORSION ANGLES, PERIOD 4\s+\(DEGREES\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_chiral_restr': r'REMARK\s+\w+\s+CHIRAL-CENTER RESTRAINTS\s+\(A\*\*3\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_gen_planes_refined': r'REMARK\s+\w+\s+GENERAL PLANES REFINED ATOMS\s+\(A\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_gen_planes_other': r'REMARK\s+\w+\s+GENERAL PLANES OTHERS\s+\(A\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_mcbond_it': r'REMARK\s+\w+\s+MAIN\-CHAIN BOND REFINED ATOMS\s+\(A\*\*2\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_mcbond_other': r'REMARK\s+\w+\s+MAIN\-CHAIN BOND OTHER ATOMS\s+\(A\*\*2\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_mcangle_it': r'REMARK\s+\w+\s+MAIN\-CHAIN ANGLE REFINED ATOMS\s+\(A\*\*2\):\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_mcangle_other': r'REMARK\s+\w+\s+MAIN\-CHAIN ANGLE OTHER ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_scbond_it': r'REMARK\s+\w+\s+SIDE\-CHAIN BOND REFINED ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_scbond_other': r'REMARK\s+\w+\s+SIDE\-CHAIN BOND OTHER ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_scangle_other': r'REMARK\s+\w+\s+SIDE\-CHAIN ANGLE OTHER ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_long_range_B_refined': r'REMARK\s+\w+\s+LONG RANGE B REFINED ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
    'r_long_range_B_other': r'REMARK\s+\w+\s+LONG RANGE B OTHER ATOMS\s+\(A\*\*2\)\s*:\s+(\d+)\s*;\s*([\.\d]+)\s*;\s*([\.\d]+)',
}

# these require matching 2 lines specified as a tuple
refine_hist_regexes = {
    'number_atoms_total': (
        r'REMARK\s+\w+\s+NUMBER OF NON\-HYDROGEN ATOMS USED IN REFINEMENT\.',
        r'REMARK\s+\w+\s+ALL ATOMS\s+:\s+(\d+)',
    )
}


def run(pdb_file, doc: cif.Document):
    refine_ls_shell_matches = OrderedDict([('pdbx_refine_id', '1')])
    refine_ls_restr_matches = OrderedDict()
    refine_hist_matches = OrderedDict([('cycle_id', '1')])
    with open(pdb_file, 'rt') as pdb:
        previous = None
        for line in pdb:
            for k, patt in refine_ls_shell_regexes.items():
                x = re.search(patt, line)
                if x:
                    refine_ls_shell_matches[k] = x.group(1)

            for k, patt in refine_ls_restr_regexes.items():
                x = re.search(patt, line)
                if x:
                    refine_ls_restr_matches[k] = (x.group(1), x.group(2), x.group(3))

            for k, patts in refine_hist_regexes.items():
                x = re.search(patts[1], line)
                if x:
                    y = re.search(patts[0], previous)
                    if y:
                        refine_hist_matches[k] = x.group(1)
            previous = line

    for k in refine_ls_shell_regexes:
        if k not in refine_ls_shell_matches:
            print(k, 'not found')
    # print('refine_ls_shell_matches:', refine_ls_shell_matches)

    for k in refine_ls_restr_regexes:
        if k not in refine_ls_restr_matches:
            print(k, 'not found')
    # print('refine_ls_restr_matches:', refine_ls_restr_matches)

    for k in refine_hist_regexes:
        if k not in refine_hist_matches:
            print(k, 'not found')
    # print('refine_hist_matches:', refine_hist_matches)

    # need to transform the values for the _refine_ls_restr loop
    type_values = []
    dev_ideal_values = []
    dev_ideal_target_values = []
    number_values = []
    for k, t in refine_ls_restr_matches.items():
        type_values.append(k)
        dev_ideal_values.append(t[1])
        dev_ideal_target_values.append(t[2])
        number_values.append(t[0])
    tx_refine_ls_restr_matches = OrderedDict(
        [
            ('pdbx_refine_id', [1] * len(type_values)),
            ('type', type_values),
            ('dev_ideal', dev_ideal_values),
            ('dev_ideal_target', dev_ideal_target_values),
            ('number', number_values),
        ]
    )

    if doc:
        pdb_deposition.create_pairs(refine_ls_shell_matches, '_' + 'refine_ls_shell.', doc[0])
        pdb_deposition.create_loop(tx_refine_ls_restr_matches, '_' + 'refine_ls_restr.', doc[0])
        pdb_deposition.create_pairs(refine_hist_matches, '_' + 'refine_hist.', doc[0])


def main():
    global LOG

    parser = argparse.ArgumentParser(description="scrape_pdb_stats")

    parser.add_argument("-m", "--mmcif", help="MMCIF to add to")
    parser.add_argument("-p", "--pdb", help="PDB to scrape")
    parser.add_argument("-o", "--output", help="MMCIF file to create")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()

    LOG = utils.Logger(logfile=args.log_file, level=args.log_level)
    utils.LOG = LOG
    LOG.info("scrape_pdb_stats: ", args)

    doc = None
    if args.pdb:
        struc = gemmi.read_pdb(args.pdb)
        doc = struc.make_mmcif_document()
    if args.mmcif:
        doc = cif.read(args.mmcif)
    else:
        doc = cif.Document()
        doc.add_new_block('xxx')

    run(args.pdb, doc)

    print(doc.as_string())


if __name__ == "__main__":
    main()
