import argparse
import re
from pathlib import Path
import gemmi
from gemmi import cif

from xchemalign import utils

from pdbdepo import pdb_deposition

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


def find_summary_section(file, begin_patt):
    count = 0
    for line in file:
        count += 1
        x = re.search(begin_patt, line)
        if x:
            return count
    return None


def match(line, pattern, keys, results):
    x = re.search(pattern, line)
    if x:
        for i, k in enumerate(keys):
            info('adding', k, i, x.group(i + 1))
            results[k] = x.group(i + 1)


def handle_summary_line(line: str, key_patt_dict: dict, reflns: dict, shell: dict):
    for k, p in key_patt_dict.items():
        if match_summary_line(line, k, p, reflns, shell):
            return


def match_summary_line(line: str, key: str, pattern: str, reflns: dict, shell: dict):
    x = re.search(pattern, line)
    if x:
        reflns[key] = x.group(1)
        # first row is outer shell, second is inner
        shell[key] = (x.group(3), x.group(2))
        return True
    else:
        return False


KEY_REFLNS_ENTRY_ID = 'entry_id'
KEY_REFLNS_DIFFRN_ID = 'pdbx_diffrn_id'
KEY_REFLNS_PDBX_ORDINAL = 'pdbx_ordinal'
KEY_REFLNS_RESO_LOW = 'd_resolution_low'
KEY_REFLNS_RESO_HIGH = 'd_resolution_high'
KEY_REFLNS_NUM_OBS = 'number_obs'
KEY_REFLNS_PDBX_RMERGE_I_OBS = 'pdbx_Rmerge_I_obs'
KEY_REFLNS_PDBX_RRIM_I_ALL = 'pdbx_Rrim_I_all'
KEY_REFLNS_PDBX_RPIM_I_ALL = 'pdbx_Rpim_I_all'
KEY_REFLNS_PDBX_CC_HALF = 'pdbx_CC_half'
KEY_REFLNS_PDBX_NUM_MEASURED = 'pdbx_number_measured_all'
KEY_REFLNS_PDBX_REDUNDANCY = 'pdbx_redundancy'
KEY_REFLNS_CHI_SQ = 'pdbx_chi_squared'
KEY_REFLNS_POSSIBLE_OBS = 'percent_possible_obs'
KEY_REFLNS_NETI_OVER_SIGMA = 'pdbx_netI_over_sigmaI'

# this is for the aimless.log files
d_autoproc = {
    KEY_REFLNS_RESO_LOW: r'Low resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_RESO_HIGH: r'High resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RMERGE_I_OBS: r'Rmerge\s+\(all\s+I\+\s*and\s*I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RRIM_I_ALL: r'Rmeas\s+\(all\s*I\+\s*&\s*I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RPIM_I_ALL: r'Rpim\s+\(all I\+\s+&\s+I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_CC_HALF: r'Mn\(I\)\s+half\-set\s+correlation\s+CC\(1/2\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_NUM_MEASURED: r'Total\s+number\s+of\s+observations\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_REDUNDANCY: r'Multiplicity\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_POSSIBLE_OBS: r'Completeness\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_NETI_OVER_SIGMA: r'Mean\(\(I\)\/sd\(I\)\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
}

# this is for the .table1 files
d_autoproc_staraniso = {
    KEY_REFLNS_RESO_LOW: r'Low resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_RESO_HIGH: r'High resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RMERGE_I_OBS: r'Rmerge\s+\(all\s+I\+\s*&\s*I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RRIM_I_ALL: r'Rmeas\s+\(all\s*I\+\s*&\s*I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RPIM_I_ALL: r'Rpim\s+\(all I\+\s+&\s+I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_CC_HALF: r'CC\(1/2\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_NUM_MEASURED: r'Total\s+number\s+of\s+observations\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_REDUNDANCY: r'Multiplicity\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_CHI_SQ: r'Mean\(Chi\^2\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_POSSIBLE_OBS: r'Completeness \(ellipsoidal\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_NETI_OVER_SIGMA: r'Mean\(I\)\/sd\(I\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
}

# this is for the .log files
d_xia_3dii = {
    KEY_REFLNS_RESO_LOW: r'Low resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_RESO_HIGH: r'High resolution limit\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RMERGE_I_OBS: r'Rmerge\s+\(all\s+I\+\s+and\s+I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RRIM_I_ALL: r'Rmeas\s+\(all\s+I\+\s+&\s+I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_RPIM_I_ALL: r'Rpim\s+\(all\s+I\+\s+&\s+I-\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_CC_HALF: r'Mn\(I\)\s+half\-set\s+correlation\s+CC\(1/2\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_NUM_MEASURED: r'Total\s+number\sof\s+observations\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_REDUNDANCY: r'Multiplicity\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_CHI_SQ: r'Mean\(Chi\^2\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_POSSIBLE_OBS: r'Completeness\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_NETI_OVER_SIGMA: r'Mean\(\(I\)\/sd\(I\)\)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
}


def handle_text_file(file, regexes):
    if file is None:
        error('Log file not defined')
        return {}, {}
    elif not Path(file).is_file():
        error('Log file ' + str(file) + ' not present')
        return {}, {}

    reflns = {KEY_REFLNS_ENTRY_ID: 'UNNAMED', KEY_REFLNS_DIFFRN_ID: 1, KEY_REFLNS_PDBX_ORDINAL: 1}
    shell = {KEY_REFLNS_DIFFRN_ID: (1, 1), KEY_REFLNS_PDBX_ORDINAL: (1, 2)}
    with open(file, "rt") as f:
        processed_count = 0
        read_count = find_summary_section(f, r'\s+Overall\s+InnerShell\s+OuterShell')
        # info('read', read_count, 'lines to find start of summary section')
        if read_count is not None:
            for line in f:
                read_count += 1
                processed_count += 1
                line = line.strip()
                if line.startswith('Anomalous'):
                    break

                handle_summary_line(line, regexes, reflns, shell)

    for k in regexes:
        if k not in reflns:
            warn('key', k, 'not found for reflns')
        if k not in shell:
            warn('key', k, 'not found for shell')

    info("read", read_count, "lines,  processed", processed_count)
    return reflns, shell


def handle_autoproc(file):
    return handle_text_file(file, d_autoproc)


def handle_autoproc_staraniso(file):
    return handle_text_file(file, d_autoproc_staraniso)


def handle_xia2_dials(file):
    return handle_text_file(file, d_xia_3dii)


def handle_xia_3dii(file):
    return handle_text_file(file, d_xia_3dii)


def handle_file(file, type, doc: cif.Document, outputfile: str):
    if type == 'autoproc':
        reflns, shell = handle_autoproc(file)
    elif type == 'xia_3dii':
        reflns, shell = handle_xia_3dii(file)
    else:
        info('Unsupported type: ' + type)
        return None

    # info('reflns:', reflns)
    # info('shell:', shell)

    if not doc:
        doc = cif.Document()
    block = doc.add_new_block('x')
    if reflns:
        pdb_deposition.create_pairs(reflns, '_reflns.', block)
    if shell:
        pdb_deposition.create_loop(shell, '_reflns_shell.', block)

    if outputfile:
        doc.write_file(outputfile)
    # else:
    #     info(doc.as_string(cif.Style.Simple))

    return doc


def main():
    global LOG

    parser = argparse.ArgumentParser(description="scrape_processing_stats")

    parser.add_argument("-f", "--file", required=True, help="Log file to parse")
    parser.add_argument("-m", "--mmcif", help="MMCIF to add to")
    parser.add_argument("-p", "--pdb", help="PDB to convert and add to")
    parser.add_argument("-o", "--output", help="MMCIF file to create")
    parser.add_argument("-t", "--type", required=True, help="Type of logfile")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()

    LOG = utils.Logger(logfile=args.log_file, level=args.log_level)
    utils.LOG = LOG
    LOG.info("scrape_processing_stats: ", args)

    doc = None
    if args.pdb:
        struc = gemmi.read_pdb(args.pdb)
        doc = struc.make_mmcif_document()
    elif args.mmcif:
        doc = cif.read(args.mmcif)

    handle_file(args.file, args.type, doc, args.output)


if __name__ == "__main__":
    main()
