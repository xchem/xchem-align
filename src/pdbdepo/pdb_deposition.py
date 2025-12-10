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
from pdbdepo import scrape_pdb_stats

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


# def gen_xtal_dir_from_pdb(pdb_file, xtal_name):
#     # not clear what is the best approach for this
#     d = Path(pdb_file)
#     while d is not None:
#         d = d.parent
#         if d.name == xtal_name:
#             return d
#     return None


def generate_dimple_log_path(xtal_dir_path):
    dimple_log_path = xtal_dir_path / 'dimple/dimple/dimple.log'
    if dimple_log_path.is_file():
        return dimple_log_path
    else:
        warn('Dimple log file not found: ' + str(dimple_log_path))
        return None


def generate_collection_info_path(xtal_dir_path, xtal_name):
    collection_info_path = xtal_dir_path / 'autoprocessing/' / (xtal_name + '_collection_info.cif')
    if collection_info_path.is_file():
        return collection_info_path
    else:
        warn('Collection info CIF file not found: ' + str(collection_info_path))
        return None


def process_input(
    base_dir: Path,
    input_path: Path,
    collator_output_dir: Path,
    soakdb_file: Path,
    meta: dict,
    mmcifgen_block,
    output_dir: Path,
):
    mmcifgen_refine_tags = None
    mmcifgen_refine_values = None
    mmcifgen_refine_item = find_loop_item(mmcifgen_block, '_refine')
    if mmcifgen_refine_item:
        mmcifgen_refine_tags = list(mmcifgen_refine_item.loop.tags)
        mmcifgen_refine_values = list(mmcifgen_refine_item.loop.values)
        mmcifgen_refine_item.erase()

    path_to_soakdb_file = base_dir / input_path / soakdb_file
    crystals = meta[Constants.META_XTALS]
    info("reading soakdb file:", path_to_soakdb_file)
    df = dbreader.read_pdb_depo(path_to_soakdb_file)
    info("read {} rows".format(len(df)))
    for index, row in df.iterrows():
        mmcif = row[Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST]
        if mmcif is None or mmcif == 'None':
            refinement_type = Constants.SOAKDB_VALUE_REFMAC
        else:
            refinement_type = Constants.SOAKDB_VALUE_BUSTER
        xtal_name = row[Constants.SOAKDB_XTAL_NAME]
        xtal_in_path = base_dir / input_path / utils.Constants.DEFAULT_MODEL_BUILDING_DIR / xtal_name
        xtal_out_path = output_dir / xtal_name
        if xtal_out_path.is_dir():
            # delete it
            shutil.rmtree(xtal_out_path)
        # create the output dir
        xtal_out_path.mkdir(parents=True)
        info(index, xtal_name, refinement_type)
        meta_data = crystals.get(xtal_name)
        if meta_data is None:
            warn("Crystal {} not found in metadata. This is strange.".format(xtal_name))
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
                        info("WARNING: no event map file for " + xtal_name)
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
                            info("Couldn't find dir ")
                            break
                        if parent.name == xtal_name:
                            mtz_free_path = parent / mtz_free
                            break
                if not (base_dir / utils.make_path_relative(mtz_free_path)).is_file():
                    warn("couldn't find mtz_free file " + mtz_free)

            info(
                "located files:\n  mtz_latest: {}\n  mtz_free: {}\n  mmcif: {}\n  pdb: {}\n  ccp4: {}".format(
                    mtz_latest, str(mtz_free_path), mmcif, pdb, ccp4_files
                )
            )

            if refinement_type == Constants.SOAKDB_VALUE_BUSTER:
                cif_doc = read_buster_structure(base_dir / utils.make_path_relative(Path(mmcif)))
            else:
                cif_doc = read_refmac_structure(base_dir / utils.make_path_relative(Path(pdb)))

            cif_block0 = cif_doc[0]
            cif_doc.write_file(str(xtal_out_path / 'original.cif'))

            # do some cleaning
            delete_pair_item(cif_doc, '_struct')
            delete_pair_item(cif_doc, '_struct_keywords')
            delete_pair_item(cif_doc, '_pdbx_database_status')

            refine_item = find_loop_item(cif_block0, '_refine')
            if refine_item:
                append_loop_items(cif_block0, refine_item, mmcifgen_refine_tags, mmcifgen_refine_values)

            # add in the common metadata (generated by mmcif-gen)
            for item in mmcifgen_block:
                cif_block0.add_item(item)

            prog = row.get(Constants.SOAKDB_COL_DATA_PROCESSING_PROGRAM)
            if prog == 'dials':
                prog = 'xia2-dials'
            if prog and prog != 'None':
                prog = prog.lower()
                info('data processing was done with ' + prog)
                logfile = row[Constants.SOAKDB_COL_DATA_PROCESSING_PATH_TO_LOGFILE]
                if logfile and logfile != 'None':
                    stats_cif = base_dir / utils.make_path_relative(Path(logfile).parent) / 'xia2.mmcif.bz2'
                    if not stats_cif.is_file():
                        p = None
                        # info(str(stats_cif) + ' mmcif file containing data processing stats not found')
                        logfile_p = Path(logfile)
                        if logfile_p.is_file():
                            p = base_dir / utils.make_path_relative(logfile_p)
                            info('processing stats from', logfile_p.name)
                        elif prog == 'autoproc' and logfile.endswith('.log'):
                            # look for the aimless.log file
                            aimless = base_dir / utils.make_path_relative(logfile_p).parent / 'aimless.log'
                            if aimless.is_file():
                                p = aimless
                                info('processing stats from aimless.log')
                            else:
                                warn('no aimless.log file found')
                        else:
                            warn('could not find log')

                        if p:
                            stats_doc = scrape_processing_stats.handle_file(p, prog, None, None)
                            if stats_doc is None:
                                warn('stats could not be scraped from log file')
                            else:
                                stats_block0 = stats_doc[0]
                                for item in stats_block0:
                                    cif_block0.add_item(item)

                    else:
                        info('processing stats from xia2.mmcif.bz2')
                        item_software = None
                        item_reflns_shell = None
                        pair_reflns = OrderedDict()
                        with bz2.open(stats_cif) as stats:
                            txt = stats.read()
                            doc = cif.read_string(txt)
                            # info("STATS:", str(doc))

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

                        existing_software_loop_item = find_loop_item(cif_block0, '_software')
                        if existing_software_loop_item is None:
                            warn('existing software loop not found')
                        elif item_software is None:
                            warn('item software loop not found')
                        else:
                            merge_loops(existing_software_loop_item.loop, item_software.loop, cif_block0)

                        if pair_reflns:
                            for k, v in pair_reflns.items():
                                cif_block0.set_pair(k, v)
                        if item_reflns_shell:
                            cif_block0.add_item(item_reflns_shell)

            # handle the collection_info.cif file
            collection_info_p = generate_collection_info_path(xtal_in_path, xtal_name)
            if collection_info_p is not None:
                doc = cif.read(str(collection_info_p))
                info('read ' + str(collection_info_p))
                for item in doc[0]:
                    cif_block0.add_item(item)

            # move coordinates to the end
            coordinates_item = find_loop_item(cif_block0, '_atom_site')
            item_count = 0
            for item in cif_block0:
                if item == coordinates_item:
                    coordinates_pos = item_count
                item_count += 1
            cif_block0.move_item(coordinates_pos, item_count - 1)

            cif_doc.write_file(str(xtal_out_path / 'structure.cif'))

            # handle the structure factors
            merge_sf.run(
                str(base_dir / utils.make_path_relative(Path(mtz_latest))),
                str(base_dir / utils.make_path_relative(mtz_free_path)),
                ccp4_files,
                str(xtal_out_path / 'structure_factors.cif'),
                output_individual=True,
            )


def delete_pair_item(doc, prefix):
    for block in doc:
        for item in block:
            if item.pair is not None:
                if item.pair[0].startswith(prefix + '.'):
                    item.erase()


def delete_loop_item(doc, prefix):
    for block in doc:
        item = find_loop_item(block, prefix)
        if item:
            print('erasing', item)
            item.erase()


def find_loop_item(block, prefix):
    for item in block:
        if item.loop is not None:
            for tag in item.loop.tags:
                if tag.startswith(prefix + '.'):
                    return item
    return None


def append_loop_items(block_to_add_to: cif.Block, item: cif.Item, tags, values):
    """Simplistic approach that merges two 1-row loops, add the data from item2 into item1 and then
    deleting item2

    :param item: Item with loop to extend
    :param tags: the tags to append
    :param values: the values to append
    :return:
    """
    # its seems that you can't just append to the existing loop's tags and value, you have to create a new loop
    old_tags = list(item.loop.tags)
    old_values = list(item.loop.values)

    item.erase()

    old_tags.extend(tags)
    old_values.extend(values)
    loop = block_to_add_to.init_loop('', old_tags)

    loop.add_row(old_values)


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
    doc = struc.make_mmcif_document()
    scrape_pdb_stats.run(str(pdb), doc)
    return doc


def create_pairs(data: dict, prefix: str, block: cif.Block):
    block.set_pairs(prefix, data)


def create_loop(data: dict, prefix: str, block: cif.Block):
    loop = block.init_loop(prefix, list(data.keys()))
    size = len(next(iter(data.values())))
    for i in range(size):
        values = []
        for k, v in data.items():
            values.append(str(v[i]))
        loop.add_row(values)


def run(collator_dir, metadata_cif, output_dir):
    meta_doc = cif.read(metadata_cif)
    meta_block = meta_doc[0]

    collator_path = Path(collator_dir)
    config = utils.read_config_file(collator_path / 'config.yaml')
    meta = utils.read_config_file(collator_path / 'meta_collator.yaml')
    base_dir = config.get(Constants.CONFIG_BASE_DIR)

    inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
    info("found {} inputs".format(len(inputs)))
    for input in inputs:
        if input[Constants.CONFIG_TYPE] == Constants.CONFIG_TYPE_MODEL_BUILDING:
            soakdb_file = utils.find_soakdb_file(input)
            input_path = input[Constants.CONFIG_DIR]
            info("Running for " + input_path)
            process_input(
                Path(base_dir), Path(input_path), collator_path, Path(soakdb_file), meta, meta_block, Path(output_dir)
            )


def main():
    # Example:
    #   python -m pdbdepo.pdb_deposition -w path_to_collator_output -o depo -l pdb_deo.log

    global LOG

    parser = argparse.ArgumentParser(description="pdb deposition")

    parser.add_argument("-w", "--collator-dir", required=True, help="collator's output dir")
    parser.add_argument("-m", "--metadata-cif", required=True, help="CIF file with common metadata (mmcif-gen)")

    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    LOG = utils.Logger(logfile=args.log_file, level=args.log_level)
    utils.LOG = LOG
    scrape_processing_stats.LOG = LOG
    merge_sf.LOG = LOG
    LOG.info("pdb_deposition: ", args)

    run(args.collator_dir, args.metadata_cif, args.output_dir)

    LOG.report()
    LOG.close()


if __name__ == "__main__":
    main()
