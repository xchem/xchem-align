# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import bz2
import shutil
from pathlib import Path
from collections import OrderedDict

import gemmi
from gemmi import cif

from xchemalign import dbreader, utils
from xchemalign.utils import Constants
from pdbdepo import merge_sf
from pdbdepo import scrape_processing_stats


def process_input(
    base_dir: Path,
    input_path: Path,
    collator_output_dir: Path,
    soakdb_file: Path,
    meta: dict,
    output_dir: Path,
    logger,
):
    path_to_soakdb_file = base_dir / input_path / soakdb_file
    crystals = meta[Constants.META_XTALS]
    logger.info("reading soakdb file:", path_to_soakdb_file)
    df = dbreader.read_pdb_depo(path_to_soakdb_file)
    logger.info("read {} rows".format(len(df)))
    for index, row in df.iterrows():
        mmcif = row[Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST]
        if mmcif is None or mmcif == 'None':
            refinement_type = Constants.SOAKDB_VALUE_REFMAC
        else:
            refinement_type = Constants.SOAKDB_VALUE_BUSTER
        xtal_name = row[Constants.SOAKDB_XTAL_NAME]
        xtal_out_dir = output_dir / xtal_name
        if xtal_out_dir.is_dir():
            # delete it
            shutil.rmtree(xtal_out_dir)
        # create the output dir
        xtal_out_dir.mkdir(parents=True)
        logger.info(index, xtal_name, refinement_type)
        meta_data = crystals.get(xtal_name)
        if meta_data is None:
            logger.warn("Crystal {} not found in metadata. This is strange.".format(xtal_name))
        else:
            mtz_latest = row.get(Constants.SOAKDB_COL_MTZ_LATEST)
            mtz_free = row.get(Constants.SOAKDB_COL_MTZ_FREE)
            pdb = row.get(Constants.SOAKDB_COL_PDB)
            crystallographic_files = meta_data.get(Constants.META_XTAL_FILES)
            ligand_binding_events = crystallographic_files.get(Constants.META_BINDING_EVENT)
            ccp4_files = []
            if ligand_binding_events:
                for ligand_binding_event in ligand_binding_events:
                    ccp4_filename = ligand_binding_event.get(Constants.META_FILE)
                    if not ccp4_filename:
                        logger.info("WARNING: no event map file for " + xtal_name)
                    else:
                        ccp4_files.append(str(collator_output_dir.parent / ligand_binding_event[Constants.META_FILE]))

            mtz_free_path = None
            if mtz_free:
                # fix the path of the mtz_free file as it is probably just the filename, not the path
                mtz_free_path = Path(mtz_free)
                if not mtz_free_path.is_absolute():
                    parent = Path(pdb)
                    while True:
                        parent = parent.parent
                        if parent == Path('/'):
                            logger.info("Couldn't find dir ")
                            break
                        if parent.name == xtal_name:
                            mtz_free_path = parent / mtz_free
                            break
                if not (base_dir / utils.make_path_relative(mtz_free_path)).is_file():
                    logger.warn("couldn't find mtz_free file " + mtz_free)

            logger.info(
                "located files:\n  mtz_latest: {}\n  mtz_free: {}\n  mmcif: {}\n  pdb: {}\n  ccp4: {}".format(
                    mtz_latest, str(mtz_free_path), mmcif, pdb, ccp4_files
                )
            )

            if refinement_type == Constants.SOAKDB_VALUE_BUSTER:
                cif_doc = read_buster_structure(base_dir / utils.make_path_relative(Path(mmcif)))
            else:
                cif_doc = read_refmac_structure(base_dir / utils.make_path_relative(Path(pdb)))

            prog = row.get(Constants.SOAKDB_COL_DATA_PROCESSING_PROGRAM)
            if prog == 'dials':
                prog = 'xia2-dials'
            if prog and prog != 'None':
                prog = prog.lower()
                logger.info('data processing was done with ' + prog)
                logfile = row[Constants.SOAKDB_COL_DATA_PROCESSING_PATH_TO_LOGFILE]
                if logfile and logfile != 'None':
                    stats_cif = base_dir / utils.make_path_relative(Path(logfile).parent) / 'xia2.mmcif.bz2'
                    if not stats_cif.is_file():
                        p = None
                        # logger.info(str(stats_cif) + ' mmcif file containing data processing stats not found')
                        logfile_p = Path(logfile)
                        if logfile_p.is_file():
                            p = base_dir / utils.make_path_relative(logfile_p)
                            logger.info('processing stats from', logfile_p.name)
                        elif prog == 'autoproc' and logfile.endswith('.log'):
                            # look for the aimless.log file
                            aimless = base_dir / utils.make_path_relative(logfile_p).parent / 'aimless.log'
                            if aimless.is_file():
                                p = aimless
                                logger.info('processing stats from aimless.log')
                            else:
                                logger.warn('no aimless.log file found')
                        else:
                            logger.warn('could not find log')

                        if p:
                            stats_doc = scrape_processing_stats.handle_file(p, prog, None, None)
                            if stats_doc is None:
                                logger.warn('stats could not be scraped from log file')
                            else:
                                cif_block0 = cif_doc[0]
                                stats_block0 = stats_doc[0]
                                for item in stats_block0:
                                    cif_block0.add_item(item)

                    else:
                        logger.info('processing stats from xia2.mmcif.bz2')
                        item_software = None
                        item_reflns_shell = None
                        pair_reflns = OrderedDict()
                        with bz2.open(stats_cif) as stats:
                            txt = stats.read()
                            doc = cif.read_string(txt)
                            # logger.info("STATS:", str(doc))

                            # grab the software loop from first block
                            for item in doc[0]:
                                if item.loop:
                                    for tag in item.loop.tags:
                                        if tag.startswith('_software'):
                                            item_software = item

                            # grab the reflns data from second block
                            for item in doc[1]:
                                if item.loop is None:
                                    p = item.pair
                                    if p[0].startswith('_reflns'):
                                        pair_reflns[p[0]] = p[1]
                                else:
                                    for tag in item.loop.tags:
                                        if tag.startswith('_reflns_shell'):
                                            item_reflns_shell = item

                        block0 = cif_doc[0]
                        existing_software_loop_item = find_loop_item(block0, '_software')
                        if existing_software_loop_item is None:
                            logger.warn('existing software loop not found')
                        elif item_software is None:
                            logger.warn('item software loop not found')
                        else:
                            merge_loops(existing_software_loop_item.loop, item_software.loop, block0)

                        if pair_reflns:
                            for k, v in pair_reflns.items():
                                block0.set_pair(k, v)
                        if item_reflns_shell:
                            block0.add_item(item_reflns_shell)

            cif_doc.write_file(str(xtal_out_dir / 'structure.mmcif'))

            # handle the structure factors
            merge_sf.run(
                str(base_dir / utils.make_path_relative(Path(mtz_latest))),
                str(base_dir / utils.make_path_relative(mtz_free_path)),
                ccp4_files,
                str(xtal_out_dir / 'structure_factors.cif'),
            )


def find_loop_item(block, starting_with):
    for item in block:
        if item.loop is not None:
            for tag in item.loop.tags:
                if tag.startswith(starting_with):
                    return item
    return None


def merge_loops(loop1: cif.Loop, loop2: cif.Loop, into: cif.Block):
    d = OrderedDict()
    for tag in loop1.tags:
        d[tag] = []
    for tag in loop2.tags:
        if tag not in d:
            d[tag] = []
    num_rows1 = int(len(loop1.values) / len(loop1.tags))
    num_values = num_rows1
    expanded_tags1 = []
    for i in range(num_rows1):
        expanded_tags1.extend(loop1.tags)
    for tag, value in zip(expanded_tags1, loop1.values):
        d[tag].append(value)
    for k, v in d.items():
        if len(v) < num_values:
            for i in range(num_rows1):
                v.append('?')

    num_rows2 = int(len(loop2.values) / len(loop2.tags))
    num_values = num_values + num_rows2
    expanded_tags2 = []
    for i in range(num_rows2):
        expanded_tags2.extend(loop2.tags)
    for tag, value in zip(expanded_tags2, loop2.values):
        d[tag].append(value)
    for k, v in d.items():
        if len(v) < num_values:
            for i in range(num_rows2):
                v.append('?')

    # fix the pdbx_ordinal column
    for k in d:
        if k.endswith('.pdbx_ordinal'):
            values = d[k]
            for i, k in enumerate(values):
                values[i] = str(i + 1)

    new_loop = into.init_loop('', list(d.keys()))
    for i in range(num_values):
        values = []
        for k, v in d.items():
            values.append(v[i])
        new_loop.add_row(values)


def prune_loop(loop, retain):
    for tag in loop.tags:
        if not tag.startswith(retain):
            loop.remove_column(tag)


def read_buster_structure(mmcif):
    doc = cif.read(str(mmcif))
    # is this a Structure or a Document?
    return doc


def read_refmac_structure(pdb):
    struc = gemmi.read_pdb(str(pdb))
    return struc.make_mmcif_document()


def run(collator_dir, output_dir, logger):
    collator_path = Path(collator_dir)
    config = utils.read_config_file(collator_path / 'config.yaml')
    meta = utils.read_config_file(collator_path / 'meta_collator.yaml')
    base_dir = config.get(Constants.CONFIG_BASE_DIR)

    inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
    logger.info("found {} inputs".format(len(inputs)))
    for input in inputs:
        if input[Constants.CONFIG_TYPE] == Constants.CONFIG_TYPE_MODEL_BUILDING:
            soakdb_file = utils.find_soakdb_file(input)
            input_path = input[Constants.CONFIG_DIR]
            logger.info("Running for " + input_path)
            process_input(
                Path(base_dir), Path(input_path), collator_path, Path(soakdb_file), meta, Path(output_dir), logger
            )


def main():
    # Example:
    #   python -m pdbdepo.pdb_deposition -w path_to_collator_output -o depo

    parser = argparse.ArgumentParser(description="pdb deposition")

    parser.add_argument("-w", "--collator-dir", required=True, help="collator's output dir")

    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    utils.LOG = logger
    scrape_processing_stats.LOG = logger
    logger.info("pdb_deposition: ", args)

    run(args.collator_dir, args.output_dir, logger)

    logger.report()
    logger.close()


if __name__ == "__main__":
    main()
