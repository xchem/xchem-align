import argparse
import re
from pathlib import Path
import gemmi
from bs4 import BeautifulSoup
from gemmi import cif

from xchemalign import utils

from pdbdepo import pdb_deposition


def info(*args, **kwargs):
    utils.log_info(*args, **kwargs)


def warn(*args, **kwargs):
    utils.log_warn(*args, **kwargs)


def error(*args, **kwargs):
    utils.log_error(*args, **kwargs)


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
    # print('handling', line)
    for k, p in key_patt_dict.items():
        if match_summary_line(line, k, p, reflns, shell):
            return


def match_summary_line(line: str, key: str, pattern: str, reflns: dict, shell: dict):
    x = re.search(pattern, line)
    if x:
        reflns[key] = x.group(1)
        # first row is outer shell, second is inner
        # hack to handle that the property name for the shell reflections are d_res_high/low
        shell[_replace_shell_key(key)] = (x.group(3), x.group(2))
        # print('    found', key,  x.group(1),  x.group(3),  x.group(2))
        return True
    else:
        return False


def _replace_shell_key(key):
    return key.replace('d_resolution_', 'd_res_').replace('number_obs', 'number_unique_obs')


KEY_REFLNS_ENTRY_ID = 'entry_id'
KEY_REFLNS_DIFFRN_ID = 'pdbx_diffrn_id'
KEY_REFLNS_PDBX_ORDINAL = 'pdbx_ordinal'
KEY_REFLNS_RESO_LOW = 'd_resolution_low'
KEY_REFLNS_RESO_HIGH = 'd_resolution_high'
KEY_REFLNS_PDBX_RMERGE_I_OBS = 'pdbx_Rmerge_I_obs'
KEY_REFLNS_PDBX_RRIM_I_ALL = 'pdbx_Rrim_I_all'
KEY_REFLNS_PDBX_RPIM_I_ALL = 'pdbx_Rpim_I_all'
KEY_REFLNS_PDBX_CC_HALF = 'pdbx_CC_half'
KEY_REFLNS_NUM_MEASURED = 'number_measured_obs'
KEY_REFLNS_NUM_OBSERVED = 'number_obs'
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
    KEY_REFLNS_NUM_MEASURED: r'Total\s+number\s+of\s+observations\s+(\d+)\s+(\d+)\s+(\d+)',
    KEY_REFLNS_NUM_OBSERVED: r'Total\s+number\s+unique\s+(\d+)\s+(\d+)\s+(\d+)',
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
    KEY_REFLNS_NUM_MEASURED: r'Total\s+number\s+of\s+observations\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_NUM_OBSERVED: r'Total\s+number\s+unique\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_PDBX_REDUNDANCY: r'Multiplicity\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
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
    KEY_REFLNS_NUM_MEASURED: r'Total\s+number\sof\s+observations\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
    KEY_REFLNS_NUM_OBSERVED: r'Total\s+number\s+unique\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)',
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
        if _replace_shell_key(k) not in shell:
            warn('key', k + '/' + _replace_shell_key(k), 'not found for shell')

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


# maps the row label in the xia2.multiplex HTML report's "Overall" table to the reflns key.
# the "Resolution" and "Completeness" rows need bespoke handling (range splitting / '%' stripping)
# so are not included here.
_XIA2_MULTIPLEX_ROW_KEYS = {
    'Observations': KEY_REFLNS_NUM_MEASURED,
    'Unique reflections': KEY_REFLNS_NUM_OBSERVED,
    'Multiplicity': KEY_REFLNS_PDBX_REDUNDANCY,
    'Mean I/σ(I)': KEY_REFLNS_NETI_OVER_SIGMA,
    'Rmerge': KEY_REFLNS_PDBX_RMERGE_I_OBS,
    'Rmeas': KEY_REFLNS_PDBX_RRIM_I_ALL,
    'Rpim': KEY_REFLNS_PDBX_RPIM_I_ALL,
    'CC½': KEY_REFLNS_PDBX_CC_HALF,
}


def _find_xia2_multiplex_overall_table(soup, panel_id='collapse_overall_All_data'):
    """Find the "Overall" stats table for the combined ("All data") dataset.

    Returns a dict of row label -> [Overall, Low resolution, High resolution] cell text,
    or None if the panel/table couldn't be found.
    """
    panel = soup.find(id=panel_id)
    if panel is None:
        return None
    table = panel.find('table')
    if table is None:
        return None
    data = {}
    for row in table.find_all('tr')[1:]:  # first row is the column header
        cells = row.find_all(['th', 'td'])
        if not cells:
            continue
        data[cells[0].get_text(strip=True)] = [c.get_text(strip=True) for c in cells[1:]]
    return data


def handle_xia2_multiplex(file):
    """Parse the xia2.multiplex HTML report (not a plain-text log, unlike the other processing
    programs), pulling stats from the "Overall" summary table of the combined "All data" dataset.
    """
    reflns = {KEY_REFLNS_ENTRY_ID: 'UNNAMED', KEY_REFLNS_DIFFRN_ID: 1, KEY_REFLNS_PDBX_ORDINAL: 1}
    shell = {KEY_REFLNS_DIFFRN_ID: (1, 1), KEY_REFLNS_PDBX_ORDINAL: (1, 2)}

    if file is None:
        error('Log file not defined')
        return reflns, shell
    elif not Path(file).is_file():
        error('Log file ' + str(file) + ' not present')
        return reflns, shell

    with open(file, 'rt', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')

    table_data = _find_xia2_multiplex_overall_table(soup)
    if table_data is None:
        warn('could not find "Overall" stats table in xia2.multiplex report ' + str(file))
        return reflns, shell

    resolution = table_data.get('Resolution (Å)')
    if resolution is None or len(resolution) != 3:
        warn('key ' + KEY_REFLNS_RESO_LOW + '/' + KEY_REFLNS_RESO_HIGH + ' not found for reflns')
    else:
        overall_low, overall_high = (v.strip() for v in resolution[0].split('-'))
        low_shell_low, low_shell_high = (v.strip() for v in resolution[1].split('-'))
        high_shell_low, high_shell_high = (v.strip() for v in resolution[2].split('-'))
        reflns[KEY_REFLNS_RESO_LOW] = overall_low
        reflns[KEY_REFLNS_RESO_HIGH] = overall_high
        # first row is outer (high resolution) shell, second is inner (low resolution) shell
        shell[_replace_shell_key(KEY_REFLNS_RESO_LOW)] = (high_shell_low, low_shell_low)
        shell[_replace_shell_key(KEY_REFLNS_RESO_HIGH)] = (high_shell_high, low_shell_high)

    for label, key in _XIA2_MULTIPLEX_ROW_KEYS.items():
        values = table_data.get(label)
        if values is None or len(values) != 3:
            warn('key ' + key + ' not found for reflns')
            continue
        reflns[key] = values[0]
        shell[_replace_shell_key(key)] = (values[2], values[1])

    completeness = table_data.get('Completeness')
    if completeness is None or len(completeness) != 3:
        warn('key ' + KEY_REFLNS_POSSIBLE_OBS + ' not found for reflns')
    else:
        reflns[KEY_REFLNS_POSSIBLE_OBS] = completeness[0].rstrip('%')
        shell[_replace_shell_key(KEY_REFLNS_POSSIBLE_OBS)] = (
            completeness[2].rstrip('%'),
            completeness[1].rstrip('%'),
        )

    return reflns, shell


def handle_file(file, type, doc: cif.Document, outputfile: str):
    if type == 'autoproc':
        reflns, shell = handle_autoproc(file)
    elif type == 'autoproc_staraniso':
        reflns, shell = handle_autoproc_staraniso(file)
    elif type == 'xia_3dii':
        reflns, shell = handle_xia_3dii(file)
    elif type == 'xia2-multiplex':
        reflns, shell = handle_xia2_multiplex(file)
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
    parser = argparse.ArgumentParser(description="scrape_processing_stats")

    parser.add_argument("-f", "--file", required=True, help="Log file to parse")
    parser.add_argument("-m", "--mmcif", help="MMCIF to add to")
    parser.add_argument("-p", "--pdb", help="PDB to convert and add to")
    parser.add_argument("-o", "--output", help="MMCIF file to create")
    parser.add_argument("-t", "--type", required=True, help="Type of logfile")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()

    LOG = utils.create_singleton_logger(args.log_file, args.log_level)
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
